module CelenoC

export predict_selenoproteins

using BioSequences
using FASTX
using GFF3
using DataFrames
using CSV
using Logging

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

# Configuration structure for external tools and resources
struct PipelineConfig
    cmsearch_exec::String
    blastp_exec::String
    tblastn_exec::String
    secis_cm_file::String
    selenoprotein_fasta::String
    selenoprotein_blast_db::String
    max_3utr_search::Int
    tblastn_evalue_cutoff::Float64
    cmsearch_evalue_cutoff::Float64
    blastp_evalue_cutoff::Float64
    cpu_cores::Int
end

# ----------------------------------------------------------------------------
# Utility Functions
# ----------------------------------------------------------------------------

"""
Checks if the required external binaries are accessible in the system's PATH.
"""
function check_executables(config::PipelineConfig)
    execs = [config.cmsearch_exec, config.blastp_exec, config.tblastn_exec]
    all_found = true
    for exec in execs
        try
            # Use `which` to find the executable
            if Sys.which(exec) === nothing
                @error "Executable not found in PATH: $exec"
                all_found = false
            end
        catch
            @error "Error checking for executable: $exec"
            all_found = false
        end
    end
    return all_found
end

"""
Checks if the required resource files exist.
"""
function check_resources(config::PipelineConfig)
    files = [config.secis_cm_file, config.selenoprotein_fasta]
    # Note: BLAST DB checking is complex; we rely on BLAST execution errors if DB is missing.
    all_found = true
    for file in files
        if !isfile(file)
            @error "Resource file not found: $file"
            all_found = false
        end
    end
    return all_found
end

# ----------------------------------------------------------------------------
# 1. Data Loading and Preprocessing
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

"""
Processes a GFF3 file to reconstruct full CDS sequences for each transcript.
This handles splicing, which is crucial for accurate translation and in-frame stop codon analysis.
"""
function reconstruct_cds_from_gff3(gff3_file::String, genome::Dict{String, LongDNA{4}})
    reconstructed_cds = Dict{String, Dict{Symbol, Any}}()
    if !isfile(gff3_file)
        @warn "GFF3 file not found. Skipping annotation-based search."
        return reconstructed_cds
    end

    reader = GFF3.Reader(open(gff3_file))
    # Temporary storage for CDS fragments grouped by parent ID
    cds_fragments = Dict{String, Vector{GFF3.Record}}()

    # 1. Collect and group CDS features
    for feature in reader
        if GFF3.featuretype(feature) == "CDS"
            # Identify the parent transcript/mRNA
            parents = GFF3.attributes(feature, "Parent")
            if isempty(parents)
                # Skip CDS features without clear parentage in this implementation
                continue
            else
                parent_id = parents[1]
            end

            if !haskey(cds_fragments, parent_id)
                cds_fragments[parent_id] = []
            end
            push!(cds_fragments[parent_id], feature)
        end
    end
    close(reader)

    # 2. Reconstruct sequences
    for (transcript_id, fragments) in cds_fragments
        if isempty(fragments) continue end

        # Check consistency
        record_id = GFF3.seqid(fragments[1])
        strand = GFF3.strand(fragments[1])

        if !haskey(genome, record_id)
            @warn "Record ID $record_id found in GFF3 but not in Genome FASTA. Skipping $transcript_id."
            continue
        end

        # Sort fragments based on biological order (5' -> 3')
        if strand == '+'
            sort!(fragments, by=GFF3.seqstart)
        elseif strand == '-'
            # For reverse strand, the biological start is the highest coordinate
             sort!(fragments, by=GFF3.seqstart, rev=true)
        else
            # Skip unstranded features
            continue
        end

        # Concatenate sequences
        full_cds = LongDNA{4}("")
        for frag in fragments
            # Basic validation during reconstruction
            if GFF3.seqid(frag) != record_id || GFF3.strand(frag) != strand
                @warn "Inconsistent strand or record ID within transcript $transcript_id. Skipping."
                full_cds = nothing
                break
            end
            seq_frag = genome[record_id][GFF3.seqrange(frag)]
            if strand == '-'
                seq_frag = reverse_complement(seq_frag)
            end
            append!(full_cds, seq_frag)
        end

        if full_cds !== nothing && length(full_cds) > 0
            # Determine the 3' boundary (genomic coordinate)
            # If +, the boundary is the highest coordinate (seqend).
            # If -, the boundary is the lowest coordinate (seqstart).
            
            # We need the min/max across all fragments, regardless of the sorting order used for reconstruction.
            if strand == '+'
                 boundary = maximum(GFF3.seqend.(fragments))
            else
                 boundary = minimum(GFF3.seqstart.(fragments))
            end

            reconstructed_cds[transcript_id] = Dict(
                :sequence => full_cds,
                :record_id => record_id,
                :strand => strand,
                :boundary => boundary
            )
        end
    end

    return reconstructed_cds
end

# ----------------------------------------------------------------------------
# 2. Annotation-Based Search
# ----------------------------------------------------------------------------

"""
Analyzes reconstructed CDS sequences to find those containing in-frame TGA codons.
"""
function find_annotation_candidates(reconstructed_cds::Dict{String, Dict{Symbol, Any}})
    candidates = Dict{String, SelenoCandidate}()

    for (transcript_id, data) in reconstructed_cds
        cds_seq = data[:sequence]
        # Check for in-frame TGA (UGA)

        # We check internal codons. We exclude the first codon (start) and the last codon (annotated stop).
        # Check from the second codon (index 4) up to the penultimate codon (length - 3).
        has_internal_tga = false
        # Ensure sequence is long enough (Start + Internal + Stop = 9bp)
        if length(cds_seq) >= 9
            for i in 4:3:length(cds_seq)-3
                if cds_seq[i:i+2] == dna"TGA"
                    has_internal_tga = true
                    break
                end
            end
        end

        if has_internal_tga
            candidate_id = "Annot_" * transcript_id
            candidate = SelenoCandidate(
                candidate_id,
                "Annotation",
                data[:record_id],
                data[:strand],
                data[:boundary],
                seq=cds_seq,
                uga=true
            )
            candidates[candidate_id] = candidate
        end
    end

    return candidates
end

# ----------------------------------------------------------------------------
# 3. Homology-Based (De Novo) Search
# ----------------------------------------------------------------------------

"""
Uses TBLASTN to find regions in the genome homologous to known selenoproteins,
and verifies the Sec (U) alignment to a Stop codon (*).
"""
function find_homology_candidates(genome_fasta::String, config::PipelineConfig)
    candidates = Dict{String, SelenoCandidate}()
    @info "Running TBLASTN (this may take some time)..."

    TBLASTN_EXEC = config.tblastn_exec
    SELENOPROTEIN_FASTA = config.selenoprotein_fasta
    TBLASTN_EVALUE_CUTOFF = config.tblastn_evalue_cutoff
    CPU_CORES = config.cpu_cores

    # Added num_threads parameter for TBLASTN
    cmd = `$TBLASTN_EXEC -query $SELENOPROTEIN_FASTA -subject $genome_fasta -evalue $TBLASTN_EVALUE_CUTOFF -num_threads $CPU_CORES -outfmt "6 std qseq sseq"`

    output = try
        read(cmd, String)
    catch e
        @error "Error executing TBLASTN: $e. Ensure BLAST+ is installed and configured correctly."
        return candidates
    end

    if isempty(output)
        @info "No TBLASTN hits found."
        return candidates
    end

    # Parse the output using CSV.jl
    blast_cols = [:qseqid, :sseqid, :pident, :length, :mismatch, :gapopen, :qstart, :qend, :sstart, :send, :evalue, :bitscore, :qseq, :sseq]
    # Use IOBuffer to treat the string output as a file
    df = try
        CSV.File(IOBuffer(output), header=blast_cols, delim='\t') |> DataFrame
    catch e
        @error "Error parsing TBLASTN output: $e"
        return candidates
    end

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
        @info "No TBLASTN hits found with the U-* alignment signature."
        return candidates
    end

    @info "Found $(nrow(df)) hits with U-* alignment signature. Clustering hits..."
    
    # Determine strand for each HSP first
    df[!, :strand] = map(df.sstart, df.send) do sstart, send
        (sstart < send) ? '+' : '-'
    end

    # Group by Query (qseqid), Subject (sseqid), and Strand to handle overlapping genes on opposite strands.
    grouped = groupby(df, [:qseqid, :sseqid, :strand])

    for group in grouped
        strand = group[1, :strand]

        # Find the most downstream coordinate across all HSPs in the group (the 3' boundary)
        if strand == '+'
            # Maximum coordinate value
            boundary = maximum(max.(group[!, :sstart], group[!, :send]))
        else
            # Minimum coordinate value
            boundary = minimum(min.(group[!, :sstart], group[!, :send]))
        end

        record_id = group[1, :sseqid]
        # Create a unique ID for this homology cluster, including strand
        candidate_id = "Homol_" * group[1, :qseqid] * "_on_" * record_id * "_" * string(boundary) * strand

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
function predict_secis(candidates::Dict{String, SelenoCandidate}, genome::Dict{String, LongDNA{4}}, config::PipelineConfig)
    if isempty(candidates) return end

    MAX_3UTR_SEARCH = config.max_3utr_search

    # 1. Extract 3' UTR sequences
    utr_sequences = Dict{String, LongDNA{4}}()

    for (id, candidate) in candidates
        if !haskey(genome, candidate.record_id) continue end
        record_seq = genome[candidate.record_id]

        # Determine 3' UTR coordinates based on strand and boundary
        if candidate.strand == '+'
            start_pos = candidate.cds_boundary + 1
            end_pos = min(length(record_seq), start_pos + MAX_3UTR_SEARCH - 1)
        elseif candidate.strand == '-'
            # Reverse strand: 3' UTR is genomically *before* the boundary
            end_pos = candidate.cds_boundary - 1
            start_pos = max(1, end_pos - MAX_3UTR_SEARCH + 1)
        else
            continue
        end

        # Check boundary conditions
        if start_pos >= end_pos || start_pos < 1 || end_pos > length(record_seq)
            continue
        end

        # Extract and orient sequence (must be 5'->3' for RNA search)
        putative_3utr = record_seq[start_pos:end_pos]
        if candidate.strand == '-'
            putative_3utr = reverse_complement(putative_3utr)
        end

        utr_sequences[id] = putative_3utr
    end

    if isempty(utr_sequences)
        @info "No valid 3' UTR sequences extracted for SECIS prediction."
        return
    end

    # 2. Use mktempdir for safe temporary file handling
    mktempdir() do tmp_dir
        tmp_fasta = joinpath(tmp_dir, "utrs.fasta")
        try
            open(FASTA.Writer, tmp_fasta) do writer
                for (id, seq) in utr_sequences
                    write(writer, FASTA.Record(id, seq))
                end
            end

            # 3. Run cmsearch on the batch
            @info "Running Infernal (cmsearch) on $(length(utr_sequences)) sequences..."

            CMSEARCH_EXEC = config.cmsearch_exec
            CPU_CORES = config.cpu_cores
            CMSEARCH_EVALUE_CUTOFF = config.cmsearch_evalue_cutoff
            SECIS_CM_FILE = config.secis_cm_file

            # Use --tblout to get structured, parsable results
            tmp_tblout = joinpath(tmp_dir, "cmsearch.tblout")
            cmd = `$CMSEARCH_EXEC --cpu $CPU_CORES --tblout $tmp_tblout -E $CMSEARCH_EVALUE_CUTOFF $SECIS_CM_FILE $tmp_fasta`

            try
                # Redirect stdout to devnull to keep the console clean, but keep stderr
                run(pipeline(cmd, stdout=devnull, stderr=stderr))
            catch e
                @error "Error executing cmsearch: $e. Ensure Infernal is installed and configured correctly."
                return
            end

            # 4. Parse the results (tblout format)
            if isfile(tmp_tblout)
                hits = readlines(tmp_tblout)
                for line in hits
                    if !startswith(line, "#")
                        # Columns are space-separated. We use regex split to handle variable whitespace in tblout
                        parts = split(line, r"\s+")
                        if length(parts) >= 1
                            candidate_id = parts[1]
                            if haskey(candidates, candidate_id)
                                # Update the candidate status
                                candidates[candidate_id].has_secis = true
                            end
                        end
                    end
                end
            end

        catch e
            @error "An error occurred during SECIS prediction: $e"
        end
        # Temp files are automatically cleaned up by mktempdir
    end
end

# ----------------------------------------------------------------------------
# 5. Confirmation and Reporting
# ----------------------------------------------------------------------------

"""
Confirms homology for annotation-based candidates using BLASTP (Efficient Batch Processing).
"""
function confirm_homology_blastp(candidates::Dict{String, SelenoCandidate}, config::PipelineConfig)
    # Filter for annotation-based candidates that passed the SECIS filter and need confirmation
    to_confirm = filter(p -> p.second.source == "Annotation" && p.second.has_secis && p.second.cds_sequence !== nothing, candidates)

    if isempty(to_confirm) return end

    @info "Running BLASTP confirmation for $(length(to_confirm)) annotation-based candidates..."

    # Use mktempdir for safer temporary file handling
    mktempdir() do tmp_dir
        tmp_fasta = joinpath(tmp_dir, "candidates.faa")
        try
            sequences_written = 0
            open(FASTA.Writer, tmp_fasta) do writer
                for (id, candidate) in to_confirm
                    # Basic validation: ensure the sequence length is a multiple of 3
                    if length(candidate.cds_sequence) % 3 != 0
                        @warn "Candidate $id CDS length is not a multiple of 3. Skipping translation."
                        continue
                    end
                    # Translate the DNA sequence (assuming standard genetic code)
                    protein_seq = try
                         translate(candidate.cds_sequence, code=standard_genetic_code)
                    catch e
                        @warn "Error during translation of candidate $id: $e. Skipping."
                        continue
                    end

                    # Crucial Step: Replace stop codons '*' (AA_Term) with 'U' (AA_U, Selenocysteine)
                    protein_seq_u = replace(protein_seq, AA_Term => AA_U)
                    
                    # Remove the final annotated stop codon if it exists (represented as U after replacement)
                    if !isempty(protein_seq_u) && protein_seq_u[end] == AA_U
                       # If the translation resulted in a terminal AA_U, it means the last codon was the annotated stop codon.
                       protein_seq_u = protein_seq_u[1:end-1]
                    end
                    
                    if !isempty(protein_seq_u)
                        write(writer, FASTA.Record(id, protein_seq_u))
                        sequences_written += 1
                    end
                end
            end

            if sequences_written == 0
                @info "No valid sequences available for BLASTP confirmation."
                return
            end

            # Run blastp on the batch
            BLASTP_EXEC = config.blastp_exec
            SELENOPROTEIN_BLAST_DB = config.selenoprotein_blast_db
            BLASTP_EVALUE_CUTOFF = config.blastp_evalue_cutoff
            CPU_CORES = config.cpu_cores

            # Added num_threads parameter for BLASTP
            cmd = `$BLASTP_EXEC -query $tmp_fasta -db $SELENOPROTEIN_BLAST_DB -evalue $BLASTP_EVALUE_CUTOFF -num_threads $CPU_CORES -outfmt 6`
            
            output = try
                read(cmd, String)
            catch e
                 @error "Error executing BLASTp: $e. Ensure BLAST+ is installed and the database ($SELENOPROTEIN_BLAST_DB) is correctly formatted."
                 return
            end

            if !isempty(output)
                # Parse results to see which IDs had hits
                lines = split(output, '\n')
                for line in lines
                    if !isempty(line)
                        parts = split(line, '\t')
                        if length(parts) >= 1
                            hit_id = parts[1]
                            if haskey(candidates, hit_id)
                                candidates[hit_id].has_homology = true
                            end
                        end
                    end
                end
            end

        catch e
            @error "An unexpected error occurred during BLASTP confirmation: $e"
        end
        # Temp files are automatically cleaned up by mktempdir
    end
end

function generate_report(candidates::Dict{String, SelenoCandidate}, config::PipelineConfig)
    @info "Generating report..."

    # Run BLASTP confirmation step
    confirm_homology_blastp(candidates, config)
    
    # Convert results to a DataFrame for easy viewing and saving
    results_df = DataFrame(
        ID = String[],
        Source = String[],
        Record_ID = String[],
        Strand = Char[],
        CDS_Boundary_Coord = Int[],
        Has_SECIS = Bool[],
        Has_Homology = Bool[],
        In_Frame_UGA = Bool[]
    )

    # Sort candidates by ID for reproducible results
    sorted_candidate_ids = sort(collect(keys(candidates)))

    for id in sorted_candidate_ids
        hit = candidates[id]
        push!(results_df, (
            hit.id,
            hit.source,
            hit.record_id,
            hit.strand,
            hit.cds_boundary,
            hit.has_secis,
            hit.has_homology,
            hit.in_frame_uga
        ))
    end

    # Filter for high confidence hits: SECIS + Homology + UGA
    high_confidence_df = filter(row -> row.Has_SECIS && row.Has_Homology && row.In_Frame_UGA, results_df)

    println("\n--- CelenoC Pipeline Report ---")
    if isempty(high_confidence_df)
        println("No high-confidence selenoprotein candidates identified.")
    else
        println("High-confidence candidates found: $(nrow(high_confidence_df))\n")
        # Display a summarized view
        show(high_confidence_df[!, [:ID, :Source, :Record_ID, :Strand, :CDS_Boundary_Coord]], allrows=true)
    end
    
    println("\nNOTE: All candidates require manual verification of the sequence alignments and gene models.")
    println("-----------------------------------\n")

    return results_df
end


# ----------------------------------------------------------------------------
# 6. Main Pipeline Orchestration
# ----------------------------------------------------------------------------

"""
    predict_selenoproteins(genome_fasta::String, annotations_gff3::String=""; kwargs...)

Runs the comprehensive selenoprotein identification pipeline.

# Arguments
- `genome_fasta`: Path to the genome sequence file (FASTA format).
- `annotations_gff3`: Optional path to the GFF3 annotation file.

# Keyword Arguments (Configuration)
- `cmsearch_exec="cmsearch"`: Path to the Infernal cmsearch executable.
- `blastp_exec="blastp"`: Path to the BLASTP executable.
- `tblastn_exec="tblastn"`: Path to the TBLASTN executable.
- `secis_cm_file`: REQUIRED. Path to the Eukaryotic SECIS Covariance Model (e.g., Rfam RF00031).
- `selenoprotein_fasta`: REQUIRED. Path to FASTA file of known selenoproteins (for TBLASTN query).
- `selenoprotein_blast_db`: REQUIRED. Path to formatted BLAST database of known selenoproteins (for BLASTP confirmation).
- `max_3utr_search=2000`: Maximum length (bp) of the 3' UTR to search for SECIS elements.
- `tblastn_evalue_cutoff=1e-5`: E-value cutoff for TBLASTN homology search.
- `cmsearch_evalue_cutoff=0.01`: E-value cutoff for cmsearch SECIS detection.
- `blastp_evalue_cutoff=1e-5`: E-value cutoff for BLASTP confirmation.
- `cpu_cores=4`: Number of CPU cores to use for external tools.
"""
function predict_selenoproteins(
    genome_fasta::String,
    annotations_gff3::String="";
    # External Tools (Defaults assume they are in PATH)
    cmsearch_exec="cmsearch",
    blastp_exec="blastp",
    tblastn_exec="tblastn",
    # Required Resources (Must be provided)
    secis_cm_file=nothing,
    selenoprotein_fasta=nothing,
    selenoprotein_blast_db=nothing,
    # Parameters
    max_3utr_search=2000,
    tblastn_evalue_cutoff=1e-5,
    cmsearch_evalue_cutoff=0.01,
    blastp_evalue_cutoff=1e-5,
    cpu_cores=4
    )

    @info "Starting CelenoC Selenoprotein Identification Pipeline..."

    # 1. Validate Configuration
    if secis_cm_file === nothing || selenoprotein_fasta === nothing || selenoprotein_blast_db === nothing
        @error "Missing required resources: secis_cm_file, selenoprotein_fasta, and selenoprotein_blast_db must be provided as keyword arguments."
        return nothing
    end

    config = PipelineConfig(
        cmsearch_exec, blastp_exec, tblastn_exec,
        secis_cm_file, selenoprotein_fasta, selenoprotein_blast_db,
        max_3utr_search, tblastn_evalue_cutoff, cmsearch_evalue_cutoff, blastp_evalue_cutoff, cpu_cores
    )

    if !check_executables(config) || !check_resources(config)
        @error "Configuration validation failed. Please check paths to executables and resource files."
        return nothing
    end

    # 2. Load Data
    if !isfile(genome_fasta)
        @error "Genome FASTA file not found: $genome_fasta"
        return nothing
    end
    
    @info "Loading genome..."
    genome = load_genome(genome_fasta)
    @info "Loaded genome with $(length(genome)) sequences."

    # 3. Initialize Candidates Dict
    all_candidates = Dict{String, SelenoCandidate}()

    # 4. Annotation-Based Search
    if !isempty(annotations_gff3) && isfile(annotations_gff3)
        @info "Starting Annotation-Based Search..."
        # Preprocess GFF3 to handle splicing
        reconstructed_cds = reconstruct_cds_from_gff3(annotations_gff3, genome)
        @info "Reconstructed $(length(reconstructed_cds)) CDS sequences from annotations."
        
        anno_candidates = find_annotation_candidates(reconstructed_cds)
        merge!(all_candidates, anno_candidates)
        @info "Found $(length(anno_candidates)) candidates with internal TGA codons."
    else
        @info "Skipping Annotation-Based Search (No valid GFF3 provided)."
    end

    # 5. Homology-Based Search
    @info "Starting Homology-Based Search..."
    homol_candidates = find_homology_candidates(genome_fasta, config)
    # Merge results
    merge!(all_candidates, homol_candidates)
    @info "Found $(length(homol_candidates)) unique candidate regions via homology search (TBLASTN)."

    if isempty(all_candidates)
        @warn "No initial candidates found (neither annotation nor homology). Exiting."
        return nothing
    end

    # 6. Unified SECIS Prediction
    @info "Starting SECIS element prediction for $(length(all_candidates)) candidates..."
    predict_secis(all_candidates, genome, config)
    secis_count = count(p -> p.second.has_secis, all_candidates)
    @info "Found SECIS elements for $secis_count candidates."

    # 7. Reporting and Final Confirmation
    results_df = generate_report(all_candidates, config)

    @info "Pipeline finished."
    return results_df
end

end # module CelenoC