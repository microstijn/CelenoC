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
        if Sys.which(exec) === nothing
            @error "Executable not found in PATH: $exec"
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

    # 1. Collect and group CDS features by Parent ID
    for feature in reader
        if GFF3.featuretype(feature) == "CDS"
            parents = GFF3.attributes(feature, "Parent")
            if isempty(parents) continue end
            
            parent_id = parents[1]
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
        record_id = GFF3.seqid(fragments[1])
        strand = GFF3.strand(fragments[1])

        if !haskey(genome, record_id) || !(strand in ['+', '-']) continue end

        # Sort fragments based on biological order (5' -> 3')
        sort!(fragments, by=GFF3.seqstart, rev=(strand == '-'))

        # Concatenate sequences
        full_cds = LongDNA{4}("")
        valid = true
        for frag in fragments
            if GFF3.seqid(frag) != record_id || GFF3.strand(frag) != strand
                valid = false; break
            end
            seq_frag = genome[record_id][GFF3.seqrange(frag)]
            if strand == '-'
                seq_frag = reverse_complement(seq_frag)
            end
            append!(full_cds, seq_frag)
        end

        if valid && length(full_cds) > 0
            # Determine the 3' boundary (genomic coordinate)
            boundary = (strand == '+') ? maximum(GFF3.seqend.(fragments)) : minimum(GFF3.seqstart.(fragments))

            reconstructed_cds[transcript_id] = Dict(:sequence => full_cds, :record_id => record_id, :strand => strand, :boundary => boundary)
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
        has_internal_tga = false
        # Check internal codons (exclude start and annotated stop).
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
            candidate = SelenoCandidate(candidate_id, "Annotation", data[:record_id], data[:strand], data[:boundary], seq=cds_seq, uga=true)
            candidates[candidate_id] = candidate
        end
    end
    return candidates
end

# ----------------------------------------------------------------------------
# 3. Homology-Based (De Novo) Search
# ----------------------------------------------------------------------------

"""
Uses TBLASTN to find regions homologous to known selenoproteins, verifying the U-* alignment.
"""
function find_homology_candidates(genome_fasta::String, config::PipelineConfig)
    candidates = Dict{String, SelenoCandidate}()
    @info "Running TBLASTN (this may take some time)..."

    # Run TBLASTN
    cmd = `$(config.tblastn_exec) -query $(config.selenoprotein_fasta) -subject $genome_fasta -evalue $(config.tblastn_evalue_cutoff) -num_threads $(config.cpu_cores) -outfmt "6 std qseq sseq"`

    output = try read(cmd, String) catch e; @error "Error executing TBLASTN: $e"; return candidates end
    if isempty(output) return candidates end

    # Parse output
    blast_cols = [:qseqid, :sseqid, :pident, :length, :mismatch, :gapopen, :qstart, :qend, :sstart, :send, :evalue, :bitscore, :qseq, :sseq]
    df = try CSV.File(IOBuffer(output), header=blast_cols, delim='\t') |> DataFrame catch e; @error "Error parsing TBLASTN output: $e"; return candidates end

    # CRITICAL: Verify U-* alignment signature
    df[!, :has_sec_signature] = map(df.qseq, df.sseq) do qseq, sseq
        any(q_char == 'U' && s_char == '*' for (q_char, s_char) in zip(qseq, sseq))
    end
    df = filter(:has_sec_signature => identity, df)

    if nrow(df) == 0 return candidates end
    @info "Found $(nrow(df)) hits with U-* signature. Clustering hits..."

    df[!, :strand] = map(df.sstart, df.send) do s, e; (s < e) ? '+' : '-' end

    # Group hits (HSPs) by Query, Subject, and Strand
    for group in groupby(df, [:qseqid, :sseqid, :strand])
        strand = group[1, :strand]
        boundary = (strand == '+') ? maximum(max.(group[!, :sstart], group[!, :send])) : minimum(min.(group[!, :sstart], group[!, :send]))
        record_id = group[1, :sseqid]
        candidate_id = "Homol_" * group[1, :qseqid] * "_on_" * record_id * "_" * string(boundary) * strand
        
        candidates[candidate_id] = SelenoCandidate(candidate_id, "Homology", record_id, strand, boundary, homology=true, uga=true)
    end
    return candidates
end

# ----------------------------------------------------------------------------
# 4. Unified SECIS Element Prediction
# ----------------------------------------------------------------------------

"""
Extracts putative 3' UTRs and runs Infernal's cmsearch efficiently on the batch.
"""
function predict_secis(candidates::Dict{String, SelenoCandidate}, genome::Dict{String, LongDNA{4}}, config::PipelineConfig)
    if isempty(candidates) return end

    utr_sequences = Dict{String, LongDNA{4}}()
    MAX_3UTR_SEARCH = config.max_3utr_search

    for (id, candidate) in candidates
        if !haskey(genome, candidate.record_id) continue end
        record_seq = genome[candidate.record_id]

        if candidate.strand == '+'
            start_pos = candidate.cds_boundary + 1
            end_pos = min(length(record_seq), start_pos + MAX_3UTR_SEARCH - 1)
        elseif candidate.strand == '-'
            end_pos = candidate.cds_boundary - 1
            start_pos = max(1, end_pos - MAX_3UTR_SEARCH + 1)
        else continue end

        if start_pos >= end_pos || start_pos < 1 || end_pos > length(record_seq) continue end

        putative_3utr = record_seq[start_pos:end_pos]
        if candidate.strand == '-'
            putative_3utr = reverse_complement(putative_3utr)
        end
        utr_sequences[id] = putative_3utr
    end

    if isempty(utr_sequences) return end

    mktempdir() do tmp_dir
        tmp_fasta = joinpath(tmp_dir, "utrs.fasta")
        open(FASTA.Writer, tmp_fasta) do writer
            for (id, seq) in utr_sequences
                write(writer, FASTA.Record(id, seq))
            end
        end

        @info "Running Infernal (cmsearch) on $(length(utr_sequences)) sequences..."
        tmp_tblout = joinpath(tmp_dir, "cmsearch.tblout")
        cmd = `$(config.cmsearch_exec) --cpu $(config.cpu_cores) --tblout $tmp_tblout -E $(config.cmsearch_evalue_cutoff) $(config.secis_cm_file) $tmp_fasta`

        try run(pipeline(cmd, stdout=devnull, stderr=stderr)) catch e; @error "Error executing cmsearch: $e"; return end

        if isfile(tmp_tblout)
            for line in readlines(tmp_tblout)
                if !startswith(line, "#")
                    parts = split(line, r"\s+")
                    if !isempty(parts) && haskey(candidates, parts[1])
                        candidates[parts[1]].has_secis = true
                    end
                end
            end
        end
    end
end

# ----------------------------------------------------------------------------
# 5. Confirmation and Reporting
# ----------------------------------------------------------------------------

"""
Confirms homology for annotation-based candidates using BLASTP.
"""
function confirm_homology_blastp(candidates::Dict{String, SelenoCandidate}, config::PipelineConfig)
    to_confirm = filter(p -> p.second.source == "Annotation" && p.second.has_secis && p.second.cds_sequence !== nothing, candidates)
    if isempty(to_confirm) return end

    @info "Running BLASTP confirmation for $(length(to_confirm)) annotation-based candidates..."
    mktempdir() do tmp_dir
        tmp_fasta = joinpath(tmp_dir, "candidates.faa")
        sequences_written = 0
        open(FASTA.Writer, tmp_fasta) do writer
            for (id, candidate) in to_confirm
                if length(candidate.cds_sequence) % 3 != 0 continue end
                
                protein_seq = try translate(candidate.cds_sequence, code=standard_genetic_code) catch e; @warn "Translation failed for $id: $e"; continue end
                
                # Crucial: Replace stop codons '*' (AA_Term) with 'U' (AA_U)
                protein_seq_u = replace(protein_seq, AA_Term => AA_U)
                
                # Remove the final annotated stop codon if present (now a U)
                if !isempty(protein_seq_u) && protein_seq_u[end] == AA_U
                   protein_seq_u = protein_seq_u[1:end-1]
                end
                
                if !isempty(protein_seq_u)
                    write(writer, FASTA.Record(id, protein_seq_u))
                    sequences_written += 1
                end
            end
        end

        if sequences_written == 0 return end

        cmd = `$(config.blastp_exec) -query $tmp_fasta -db $(config.selenoprotein_blast_db) -evalue $(config.blastp_evalue_cutoff) -num_threads $(config.cpu_cores) -outfmt 6`
        
        output = try read(cmd, String) catch e; @error "Error executing BLASTp: $e"; return end

        if !isempty(output)
            for line in split(output, '\n')
                if !isempty(line)
                    parts = split(line, '\t')
                    if !isempty(parts) && haskey(candidates, parts[1])
                        candidates[parts[1]].has_homology = true
                    end
                end
            end
        end
    end
end

function generate_report(candidates::Dict{String, SelenoCandidate}, config::PipelineConfig)
    @info "Generating report..."
    confirm_homology_blastp(candidates, config)
    
    results_df = DataFrame(ID = String[], Source = String[], Record_ID = String[], Strand = Char[], CDS_Boundary_Coord = Int[], Has_SECIS = Bool[], Has_Homology = Bool[], In_Frame_UGA = Bool[])

    for id in sort(collect(keys(candidates)))
        hit = candidates[id]
        push!(results_df, (hit.id, hit.source, hit.record_id, hit.strand, hit.cds_boundary, hit.has_secis, hit.has_homology, hit.in_frame_uga))
    end

    high_confidence_df = filter(row -> row.Has_SECIS && row.Has_Homology && row.In_Frame_UGA, results_df)

    println("\n--- CelenoC Prediction Report ---")
    if isempty(high_confidence_df)
        println("No high-confidence selenoprotein candidates identified.")
    else
        println("High-confidence candidates found: $(nrow(high_confidence_df))\n")
        show(high_confidence_df[!, [:ID, :Source, :Record_ID, :Strand, :CDS_Boundary_Coord]], allrows=true)
    end
    println("\nNOTE: Manual verification of alignments and gene models is required.")
    println("-----------------------------------\n")
    return results_df
end

# ----------------------------------------------------------------------------
# 6. Main Pipeline Orchestration
# ----------------------------------------------------------------------------

"""
    predict_selenoproteins(genome_fasta::String, annotations_gff3::String=""; kwargs...)

Runs the comprehensive de novo selenoprotein identification pipeline.

This pipeline identifies candidates based on three key pieces of evidence:
1.  An in-frame UGA (TGA) codon, found either via annotation or homology.
2.  A downstream SECIS RNA element in the 3' UTR.
3.  Homology to a known selenoprotein family.

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
    # External Tools
    cmsearch_exec="cmsearch",
    blastp_exec="blastp",
    tblastn_exec="tblastn",
    # Required Resources
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

    @info "Starting CelenoC Selenoprotein Prediction Pipeline..."

    # 1. Validation
    if secis_cm_file === nothing || selenoprotein_fasta === nothing || selenoprotein_blast_db === nothing
        @error "Missing required resources: secis_cm_file, selenoprotein_fasta, and selenoprotein_blast_db must be provided."
        return nothing
    end

    config = PipelineConfig(
        cmsearch_exec, blastp_exec, tblastn_exec,
        secis_cm_file, selenoprotein_fasta, selenoprotein_blast_db,
        max_3utr_search, tblastn_evalue_cutoff, cmsearch_evalue_cutoff, blastp_evalue_cutoff, cpu_cores
    )

    if !check_executables(config) || !check_resources(config)
        @error "Configuration validation failed. Check paths to executables and resource files."
        return nothing
    end

    # 2. Load Data
    if !isfile(genome_fasta)
        @error "Genome FASTA file not found: $genome_fasta"
        return nothing
    end
    @info "Loading genome..."
    genome = load_genome(genome_fasta)

    # 3. Initialize Candidates
    all_candidates = Dict{String, SelenoCandidate}()

    # 4. Annotation-Based Search
    if !isempty(annotations_gff3) && isfile(annotations_gff3)
        @info "Starting Annotation-Based Search..."
        reconstructed_cds = reconstruct_cds_from_gff3(annotations_gff3, genome)
        anno_candidates = find_annotation_candidates(reconstructed_cds)
        merge!(all_candidates, anno_candidates)
        @info "Found $(length(anno_candidates)) candidates via annotation analysis."
    else
        @info "Skipping Annotation-Based Search (No GFF3 provided)."
    end

    # 5. Homology-Based Search
    @info "Starting Homology-Based Search..."
    homol_candidates = find_homology_candidates(genome_fasta, config)
    merge!(all_candidates, homol_candidates)
    @info "Found $(length(homol_candidates)) candidates via homology search."

    if isempty(all_candidates)
        @warn "No initial candidates found. Exiting."
        return nothing
    end

    # 6. Unified SECIS Prediction
    @info "Starting SECIS element prediction for $(length(all_candidates)) candidates..."
    predict_secis(all_candidates, genome, config)

    # 7. Reporting and Final Confirmation
    results_df = generate_report(all_candidates, config)
    @info "Pipeline finished."
    return results_df
end

end # module CelenoC