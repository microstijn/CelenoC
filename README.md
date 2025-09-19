# Selenoprotein identification pipeline

## not currently maintained nor responding to issues. 

A hybrid pipeline implemented in Julia for identifying eukaryotic selenocysteine-containing proteins (selenoproteins) in genomic datasets. It combines annotation-based analysis with *de novo* homology-based searching.

## Table of Contents
1. [Introduction](#introduction)
2. [Pipeline Methodology](#pipeline-methodology)
3. [Prerequisites](#prerequisites)
4. [Installation of External Tools](#installation-of-external-tools)
5. [Configuration and Data Requirements](#configuration-and-data-requirements)
6. [Usage](#usage)
7. [Limitations](#limitations)

## Introduction
Selenoproteins incorporate the amino acid selenocysteine (Sec, U) via a complex mechanism that repurposes the UGA stop codon. This often leads to misannotations or completely missed genes. This pipeline integrates sequence analysis, RNA structure prediction, and comprehensive homology searching to identify known and novel putative selenoprotein genes.

## Pipeline Methodology
### Phase 1: Candidate Identification
The pipeline identifies candidate regions using two parallel approaches:

1.  **Annotation-Based Search (Optional):** If a GFF3 file is provided, the pipeline scans existing CDS features for in-frame UGA codons.
2.  **Homology-Based Search (De Novo):** The pipeline runs **TBLASTN** to align a database of known selenoproteins against the target genome. Crucially, it analyzes the alignments to identify the hallmark signature of selenoproteins: a Sec (U) in the query protein aligning to a Stop codon (*) in the translated genome.

### Phase 2: Unified SECIS Prediction (Batch Mode)
Candidates from both approaches are merged. The pipeline determines the putative 3' Untranslated Region (3' UTR) downstream of the annotation boundary or the homology hit boundary.

These 3' UTR sequences are batched together and analyzed efficiently in a single execution of **Infernal** (`cmsearch`) using a Covariance Model (CM) of the eukaryotic SECIS element.

### Phase 3: Integration and Confirmation
The results are integrated.

*   Homology-based hits (that possess the U-* signature) with a downstream SECIS element are reported as potential novel selenoproteins.
*   Annotation-based hits with an in-frame UGA and a SECIS element proceed to confirmation. They are translated (UGA->U) and analyzed via **BLASTP** (in batch mode) against known selenoproteins.

## Prerequisites

### Julia Environment
*   **Julia:** Version 1.6 or later. (Download from [https://julialang.org/downloads/](https://julialang.org/downloads/))
*   **Julia Packages:**
    *   `BioSequences`
    *   `FASTX`
    *   `GFF3`
    *   `DataFrames`
    *   `CSV` (For parsing BLAST results)

You can install the Julia packages using the Julia REPL:
```julia
using Pkg
Pkg.add(["BioSequences", "FASTX", "GFF3", "DataFrames", "CSV"])
```

### External Software Dependencies
The pipeline relies on external bioinformatics tools, which must be installed locally and accessible in your system's PATH.

1.  **Infernal (Infrastructure for RNA Localization)**
	    *   Required tool: `cmsearch`
	    *   Purpose: SECIS element detection.

2.  **NCBI BLAST+**
	    *   Required tools: `blastp`, `tblastn` (and `makeblastdb` for setup)
	    *   Purpose: `tblastn` is used for *de novo* homology search; `blastp` is used for confirmation of annotation-based candidates.

## Installation of External Tools

If you are using the Conda package manager, installation is straightforward.

### Infernal
*   **Linux/macOS (via Conda):**
```bash
conda install -c bioconda infernal
```
*   **Manual Installation:** Download from the [Eddy Lab website](http://eddylab.org/infernal/).

Ensure `cmsearch` is accessible:
```bash
cmsearch -h
```

### NCBI BLAST+
*   **Linux/macOS (via Conda):**
```bash
conda install -c bioconda blast
```
*   **Manual Installation:** Download from the [NCBI FTP site](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).

Ensure `blastp` and `tblastn` are accessible:
```bash
blastp -h
tblastn -h
```

## Configuration and Data Requirements
Before running the pipeline, you must prepare the required data files and configure the paths within the Julia script (e.g., `pipeline.jl`).

### 1. SECIS Covariance Model (CM)
*   **Download:** Obtain the eukaryotic model ([RF00031](https://rfam.org/family/RF00031)) from the Rfam database.

### 2. Selenoprotein Databases
You need a curated set of known selenoproteins. **Crucially, ensure Selenocysteine is represented by the letter 'U' in these files.**

*   **FASTA format (for TBLASTN query):** A standard protein FASTA file (e.g., `selenoproteins.fasta`).
*   **BLAST DB format (for BLASTP confirmation):** Format the FASTA file using `makeblastdb`:
```bash
	makeblastdb -in selenoproteins.fasta -dbtype prot -out selenodb_blast_db
```

### 3. Update Paths in the Script
Open the Julia script and modify the configuration section:

```julia
# ----------------------------------------------------------------------------
# Configuration: Paths to external tools and resources
# ----------------------------------------------------------------------------

# If installed in PATH, these defaults should work:
const CMSEARCH_EXEC = "cmsearch"
const BLASTP_EXEC = "blastp"
const TBLASTN_EXEC = "tblastn"
	
# Configuration paths (YOU MUST UPDATE THESE PATHS)
const SECIS_CM_FILE = "/path/to/Eukaryotic_SECIS.cm"
	
# Known selenoproteins (FASTA format for TBLASTN query)
const SELENOPROTEIN_FASTA = "/path/to/selenoproteins.fasta"
	
# Formatted BLAST DB for BLASTP confirmation
const SELENOPROTEIN_BLAST_DB = "/path/to/selenodb_blast_db"
	
# Parameters
const MAX_3UTR_SEARCH = 2000 
const TBLASTN_EVALUE_CUTOFF = 1e-5
const CMSEARCH_EVALUE_CUTOFF = 0.01
const CPU_CORES = 4 # Adjust based on available resources
```

## Usage
The pipeline requires the target organism's genome sequence (FASTA) and optionally the gene annotations (GFF3).

1.  **Genome Sequence:** A FASTA file (`.fna`, `.fasta`).
2.  **Gene Annotations (Optional):** A GFF3 file (`.gff`, `.gff3`).

Load the script into the Julia REPL and execute it:

```julia
# Load the script
include("pipeline.jl")
	
genome_file = "/path/to/your_genome.fna"
annotation_file = "/path/to/your_annotations.gff3"
	
# Run comprehensive pipeline (Annotation + Homology)
# If annotation_file path is invalid or empty, it will skip the annotation-based search.
run_pipeline(genome_file, annotation_file)
	
# Run only de novo pipeline (Homology only)
# run_pipeline(genome_file)
```

## Limitations
*   **Manual Verification is Essential:** This pipeline identifies high-confidence candidates, but manual inspection of the alignments (both TBLASTN U-* alignments and BLASTP U-U/C alignments) is mandatory for final confirmation.
*   **Gene Model Accuracy (De Novo):** The homology-based approach identifies regions of interest and defines the 3' boundary based on TBLASTN hits. It does not perform full gene prediction (e.g., identifying exact splice sites or the start codon). The resulting hits require manual curation or further analysis with tools like Exonerate or GeneWise to finalize the gene model.
*   **Splicing (Annotation-Based):** The annotation-based search uses a simplified GFF3 parser that may not correctly reconstruct the full CDS sequence for complex, multi-exon (spliced) genes if they are analyzed feature-by-feature.
*   **Eukaryotic Focus:** The search strategy (3' UTR location and CM model) is tailored for eukaryotes.
*   **Resource Intensity:** TBLASTN and `cmsearch` can be computationally intensive on large genomes.
