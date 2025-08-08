
# ----------------------------------------------------------------------------
# preamble
# ----------------------------------------------------------------------------

using BioSequences
using FASTX
using GFF3
using DataFrames
using CSV

# ----------------------------------------------------------------------------
# Configuration: Paths to external tools and resources
# ----------------------------------------------------------------------------

# Ensure Infernal (cmsearch) and BLAST+ (blastp, tblastn) are installed and accessible
const CMSEARCH_EXEC = "cmsearch"
const BLASTP_EXEC = "blastp"
const TBLASTN_EXEC = "tblastn"

# Configuration paths (Update these paths)
# Covariance Model for Eukaryotic SECIS (e.g., from Rfam RF00031)
const SECIS_CM_FILE = "/path/to/Eukaryotic_SECIS.cm"

# BLAST database of known selenoproteins (FASTA format for TBLASTN query)
# Ensure 'U' is used for Selenocysteine in this file.
const SELENOPROTEIN_FASTA = "/path/to/selenoproteins.fasta"
# BLAST database for confirmation (Formatted BLAST DB for BLASTP)
const SELENOPROTEIN_BLAST_DB = "/path/to/selenodb_blast_db"


const MAX_3UTR_SEARCH = 2000 # Max length of 3' UTR to search
const TBLASTN_EVALUE_CUTOFF = 1e-5
const CMSEARCH_EVALUE_CUTOFF = 0.01
const CPU_CORES = 4 # Adjust based on available resources

# ----------------------------------------------------------------------------
# Data Structures
# ----------------------------------------------------------------------------

# Represents a potential selenoprotein gene location
 mutable struct SelenoCandidate
    id::String
    source::String # "Annotation" or "Homology"
    record_id::String
    strand::Char
    # Coordinates defining the 3' boundary of the CDS/Homology hit
    cds_boundary::Int
    # Store sequence for annotation-based candidates for validation
    cds_sequence::Union{LongDNA{4}, Nothing}
    has_secis::Bool
    has_homology::Bool
    in_frame_uga::Bool # True if derived from annotation check or TBLASTN U-* alignment

    SelenoCandidate(id, source, record_id, strand, boundary; seq=nothing, uga=false, homology=false) =
        new(id, source, record_id, strand, boundary, seq, false, homology, uga)
end

# ----------------------------------------------------------------------------
# 1. Data Loading
# ----------------------------------------------------------------------------

function load_genome(fasta_file::String)
    genome = Dict{String, LongDNA{4}}()
    reader = FASTA.Reader(open(fasta_file))
    for record in reader
        genome[FASTA.identifier(record)] = FASTA.sequence(LongDNA{4}, record)
    end
    close(reader)
    return genome
end

# ----------------------------------------------------------------------------
# 2. Annotation-Based Search
# ----------------------------------------------------------------------------

"""
Analyzes GFF3 annotations (if provided) to find CDS features containing in-frame TGA codons.
"""
function find_annotation_candidates(gff3_file::String, genome::Dict{String, LongDNA{4}})
    candidates = Dict{String, SelenoCandidate}()
    if !isfile(gff3_file)
        println("GFF3 file not found. Skipping annotation-based search.")
        return candidates
    end

    reader = GFF3.Reader(open(gff3_file))

    # Note: GFF3 parsing can be complex. This implementation simplifies by analyzing individual CDS features.
    # A production pipeline should ideally reconstruct transcripts from exons.
    for feature in reader
        if GFF3.featuretype(feature) == "CDS"
            record_id = GFF3.seqid(feature)
            if !haskey(genome, record_id) continue end

            # Extract sequence
            cds_seq = genome[record_id][GFF3.seqrange(feature)]
            strand = GFF3.strand(feature)
            if strand == '-'
                cds_seq = reverse_complement(cds_seq)
            end

            # Check for in-frame TGA (UGA)
            has_potential = false
            # Check internal codons (excluding the final annotated stop codon if present in the feature)
            # We check up to length-5 to ensure we aren't checking the very last codon.
            for i in 1:3:length(cds_seq)-5
                if cds_seq[i:i+2] == dna"TGA"
                    has_potential = true
                    break
                end
            end

            if has_potential
                # Define the 3' boundary
                boundary = (strand == '+') ? GFF3.seqend(feature) : GFF3.seqstart(feature)
                # Generate a unique ID for this CDS feature
                # Attempt to retrieve a meaningful name/ID
                attrs = GFF3.attributes(feature)
                name = get(attrs, "Name", get(attrs, "ID", "CDS_" * string(GFF3.seqstart(feature))))
                candidate_id = "Annot_" * record_id * "_" * name

                # Avoid overwriting if multiple CDS segments belong to the same gene (simplification)
                if !haskey(candidates, candidate_id)
                    candidate = SelenoCandidate(candidate_id, "Annotation", record_id, strand, boundary, seq=cds_seq, uga=true)
                    candidates[candidate_id] = candidate
                end
            end
        end
    end
    close(reader)
    return candidates
end

# ----------------------------------------------------------------------------
# 3. Homology-Based (De Novo) Search
# ----------------------------------------------------------------------------

"""
Uses TBLASTN to find regions in the genome homologous to known selenoproteins, 
and verifies the Sec (U) alignment to a Stop codon (*).
"""
function find_homology_candidates(genome_fasta::String)
    candidates = Dict{String, SelenoCandidate}()
    println("Running TBLASTN (this may take some time)...")

    # Run TBLASTN: Query=Known Selenoproteins, Subject=Target Genome
    # We request standard output + qseq (query sequence) and sseq (subject/translated sequence)
    # Using -subject is simpler than creating a DB, suitable for a single run.
    cmd = `$TBLASTN_EXEC -query $SELENOPROTEIN_FASTA -subject $genome_fasta -evalue $TBLASTN_EVALUE_CUTOFF -outfmt "6 std qseq sseq"`
    
    output = try
        read(cmd, String)
    catch e
        @error "Error executing TBLASTN: $e"
        return candidates
    end

    if isempty(output) return candidates end

    # Parse the output using CSV.jl
    blast_cols = [:qseqid, :sseqid, :pident, :length, :mismatch, :gapopen, :qstart, :qend, :sstart, :send, :evalue, :bitscore, :qseq, :sseq]
    # Use IOBuffer to treat the string output as a file
    df = CSV.File(IOBuffer(output), header=blast_cols, delim='\t') |> DataFrame

    # CRITICAL STEP: Verify the U-* alignment signature
    df[!, :has_sec_signature] = falses(nrow(df))

    for i in 1:nrow(df)
        qseq = df[i, :qseq]
        sseq = df[i, :sseq]
        # Iterate through the aligned sequences
        for (q_char, s_char) in zip(qseq, sseq)
            # Check if Selenocysteine (U) in the query aligns with a Stop codon (*) in the subject
            if q_char == 'U' && s_char == '*'
                df[i, :has_sec_signature] = true
                break
            end
        end
    end

    # Filter out hits that lack the signature (e.g., Cys-homologs that don't align U to *)
    df = filter(:has_sec_signature => identity, df)
    
    if nrow(df) == 0
        println("No TBLASTN hits found with the U-* alignment signature.")
        return candidates
    end

    println("Found $(nrow(df)) hits with U-* alignment signature. Clustering hits...")

    # Process hits to define the search area for the 3' UTR.
    # Group by Query (qseqid) and Subject (sseqid) to handle fragmented hits (HSPs).
    grouped = groupby(df, [:qseqid, :sseqid])

    for group in grouped
        # Determine strand based on coordinates (sstart vs send) of the first HSP
        sstart1 = group[1, :sstart]
        send1 = group[1, :send]
        strand = (sstart1 < send1) ? '+' : '-'

        # Find the most downstream coordinate across all HSPs in the group (the 3' boundary)
        if strand == '+'
            # Maximum coordinate value
            boundary = maximum(max.(group[!, :sstart], group[!, :send]))
        else
            # Minimum coordinate value
            boundary = minimum(min.(group[!, :sstart], group[!, :send]))
        end

        record_id = group[1, :sseqid]
        # Create a unique ID for this homology cluster
        candidate_id = "Homol_" * group[1, :qseqid] * "_on_" * record_id * "_" * string(boundary)

        # We don't extract the CDS sequence here as the exact gene model is unknown
        # Homology is true by definition, UGA is true because we checked the U-* signature.
        candidate = SelenoCandidate(candidate_id, "Homology", record_id, strand, boundary, homology=true, uga=true)
        candidates[candidate_id] = candidate
    end

    return candidates
end

# ----------------------------------------------------------------------------
# 4. Unified SECIS Element Prediction (Efficient Batch Processing)
# ----------------------------------------------------------------------------

"""
Extracts putative 3' UTRs for all candidates, writes them to a single FASTA file, 
and runs Infernal's cmsearch efficiently on the batch.
"""
function predict_secis(candidates::Dict{String, SelenoCandidate}, genome::Dict{String, LongDNA{4}})
    if isempty(candidates) return end

    # 1. Extract 3' UTR sequences
    utr_sequences = Dict{String, LongDNA{4}}()

    for (id, candidate) in candidates
        if !haskey(genome, candidate.record_id) continue end
        record_seq = genome[candidate.record_id]

        # Determine 3' UTR coordinates based on strand and boundary (Annotation or Homology boundary)
        if candidate.strand == '+'
            start_pos = candidate.cds_boundary + 1
            end_pos = min(length(record_seq), start_pos + MAX_3UTR_SEARCH - 1)
        else
            # Reverse strand: 3' UTR is genomically *before* the boundary
            end_pos = candidate.cds_boundary - 1
            start_pos = max(1, end_pos - MAX_3UTR_SEARCH + 1)
        end

        if start_pos >= end_pos || start_pos < 1 || end_pos > length(record_seq) continue end

        # Extract and orient sequence (must be 5'->3' for RNA search)
        putative_3utr = record_seq[start_pos:end_pos]
        if candidate.strand == '-'
            putative_3utr = reverse_complement(putative_3utr)
        end
        
        utr_sequences[id] = putative_3utr
    end

    if isempty(utr_sequences) return end

    # 2. Write UTRs to a temporary FASTA file
    tmp_fasta = tempname() * ".fasta"
    try
        open(FASTA.Writer, tmp_fasta) do writer
            for (id, seq) in utr_sequences
                write(writer, FASTA.Record(id, seq))
            end
        end

        # 3. Run cmsearch on the batch
        println("Running Infernal (cmsearch) on $(length(utr_sequences)) sequences...")
        # Use --tblout to get structured, parsable results
        tmp_tblout = tempname() * ".tblout"
        # Use -E to apply the E-value threshold during the search for efficiency
        cmd = `$CMSEARCH_EXEC --cpu $CPU_CORES --tblout $tmp_tblout -E $CMSEARCH_EVALUE_CUTOFF $SECIS_CM_FILE $tmp_fasta`
        
        try
            # Redirect stdout to devnull to keep the console clean, but keep stderr
            run(pipeline(cmd, stdout=devnull, stderr=stderr))
        catch e
            @error "Error executing cmsearch: $e"
            return
        end

        # 4. Parse the results (tblout format)
        if isfile(tmp_tblout)
            hits = readlines(tmp_tblout)
            for line in hits
                if !startswith(line, "#")
                    # Columns are space-separated. We need the target name (col 1)
                    parts = split(line)
                    if length(parts) >= 1
                        candidate_id = parts[1]
                        if haskey(candidates, candidate_id)
                            # Update the candidate status
                            candidates[candidate_id].has_secis = true
                        end
                    end
                end
            end
            rm(tmp_tblout)
        end

    finally
        # Clean up temporary FASTA file
        isfile(tmp_fasta) && rm(tmp_fasta)
    end
end

# ----------------------------------------------------------------------------
# 5. Confirmation and Reporting
# ----------------------------------------------------------------------------

"""
Confirms homology for annotation-based candidates using BLASTP (Efficient Batch Processing).
(Homology-based candidates inherently have homology and do not need this step).
"""
function confirm_homology_blastp(candidates::Dict{String, SelenoCandidate})
    # Filter for annotation-based candidates that passed the SECIS filter and need confirmation
    to_confirm = filter(p -> p.second.source == "Annotation" && p.second.has_secis && p.second.cds_sequence !== nothing, candidates)
    
    if isempty(to_confirm) return end

    println("Running BLASTP confirmation for $(length(to_confirm)) annotation-based candidates...")

    # Write protein sequences to temp file
    tmp_fasta = tempname() * ".fasta"
    try
        open(FASTA.Writer, tmp_fasta) do writer
            for (id, candidate) in to_confirm
                 # Basic validation: ensure the sequence length is a multiple of 3 before translation
                if length(candidate.cds_sequence) % 3 != 0
                    # This might happen due to simplified GFF parsing of spliced genes
                    continue
                end
                # Translate the DNA sequence
                protein_seq = translate(candidate.cds_sequence, code=standard_genetic_code)
                # Crucial Step: Replace stop codons '*' (AA_Term) with 'U' (AA_U, Selenocysteine)
                protein_seq_u = replace(protein_seq, AA_Term => AA_U)
                write(writer, FASTA.Record(id, protein_seq_u))
            end
        end

        # Run blastp on the batch (using tabular output format 6)
        cmd = `$BLASTP_EXEC -query $tmp_fasta -db $SELENOPROTEIN_BLAST_DB -evalue 1e-5 -outfmt 6`
        output = read(cmd, String)

        if !isempty(output)
            # Parse results to see which IDs had hits
            lines = split(output, '\n')
            for line in lines
                if !isempty(line)
                    hit_id = split(line, '\t')[1]
                    if haskey(candidates, hit_id)
                        candidates[hit_id].has_homology = true
                        # Note: Manual alignment verification (U aligning to U/C) is still highly recommended.
                    end
                end
            end
        end

    catch e
        @error "Error executing BLASTp: $e"
    finally
        isfile(tmp_fasta) && rm(tmp_fasta)
    end
end

function generate_report(candidates::Dict{String, SelenoCandidate})
    println("\n--- Pipeline Report ---")
    
    # Run BLASTP confirmation step specifically for annotation hits (if any)
    confirm_homology_blastp(candidates)

    # Filter for high-confidence hits
    # Criteria: Must have SECIS AND Homology AND evidence of in-frame UGA (checked during candidate generation)
    final_hits = filter(p -> p.second.has_secis && p.second.has_homology && p.second.in_frame_uga, candidates)

    if isempty(final_hits)
        println("No high-confidence selenoprotein candidates identified.")
        return
    end

    # Separate by source for reporting clarity
    anno_hits = [c.second for (id, c) in final_hits if c.second.source == "Annotation"]
    homol_hits = [c.second for (id, c) in final_hits if c.second.source == "Homology"]

    println("\nAnnotation-Based Candidates (Confirmation of existing annotations):")
    if isempty(anno_hits)
        println("- None")
    else
        for hit in anno_hits
            println("- ID: $(hit.id), Location: $(hit.record_id) ($(hit.strand)) near $(hit.cds_boundary)")
        end
    end

    println("\nHomology-Based Candidates (Potential novel selenoproteins):")
    if isempty(homol_hits)
        println("- None")
    else
        for hit in homol_hits
             println("- ID: $(hit.id), Location: $(hit.record_id) ($(hit.strand)) near $(hit.cds_boundary)")
        end
    end

    println("\nTotal high-confidence candidates: $(length(final_hits))")
    println("NOTE: All candidates require manual verification of the sequence alignments.")
end


# ----------------------------------------------------------------------------
# 6. Main Pipeline Orchestration
# ----------------------------------------------------------------------------

"""
Main function to run the comprehensive pipeline.
genome_fasta: Required path to the genome sequence file.
annotations_gff3: Optional path to the GFF3 annotation file.
"""
function run_pipeline(genome_fasta::String, annotations_gff3::String="")
    println("Starting Comprehensive Selenoprotein Identification Pipeline...")

    # 1. Load Data
    if !isfile(genome_fasta)
        @error "Genome FASTA file not found: $genome_fasta"
        return
    end
    genome = load_genome(genome_fasta)

    # 2. Initialize Candidates Dict
    all_candidates = Dict{String, SelenoCandidate}()

    # 3. Annotation-Based Search
    if !isempty(annotations_gff3)
        anno_candidates = find_annotation_candidates(annotations_gff3, genome)
        merge!(all_candidates, anno_candidates)
        println("Found $(length(anno_candidates)) candidates via annotation analysis.")
    else
        println("Skipping Annotation-Based Search (No GFF3 provided).")
    end

    # 4. Homology-Based Search
    homol_candidates = find_homology_candidates(genome_fasta)
    # Merge results into the main dictionary
    merge!(all_candidates, homol_candidates)
    println("Found $(length(homol_candidates)) unique candidate regions via homology search (TBLASTN).")

    if isempty(all_candidates)
        println("No initial candidates found. Exiting.")
        return
    end

    # 5. Unified SECIS Prediction
    predict_secis(all_candidates, genome)

    # 6. Reporting and Final Confirmation
    generate_report(all_candidates)

    println("Pipeline finished.")
end

# Example usage (requires input files and external tools configured):
# To run both modes (if GFF3 is available):
# run_pipeline("genome.fna", "annotations.gff3") 
# To run only de novo (homology) mode:
# run_pipeline("genome.fna")