This is a selection of scripts I use to trawl through the NCBI SRA database using a reference sequence to find matches in sequencing reads that may not have made it into the respective assemblies.

# Scripts Overview

This repository contains three main scripts that work together to identify and extract plasmid sequences from sequencing data:

1. **`sra-trawler.sh`** - Downloads and screens SRA entries for potential plasmid (or other) sequences
2. **`plasmid_extractor.sh`** - Iteratively extracts and refines plasmid sequences from sequencing reads
3. **`plasmid_optimizer.sh`** - Optimizes extracted plasmids by removing repeats and adjusting coordinates

## 1. SRA Trawler (`sra-trawler.sh`)

The SRA Trawler script is designed to systematically search through the NCBI SRA database to find sequencing datasets that contain potential plasmid sequences matching a supplied reference.

### Purpose
- Downloads all SRA entries for a given organism queries or CSV input
- Maps downloaded sequencing reads to a reference genome
- Identifies entries with sufficient coverage (>1x) as potential plasmid sources
- Saves only reads that map to the reference for downstream analysis
- Maintains a SQLite database to track processing status


### Usage Examples
```bash
# Create new database from SRA query (default: fungi)
./sra-trawler.sh -f reference.fasta

# Create new database for specific organism
./sra-trawler.sh -o "Aspergillus" -f reference.fasta

# Create new database from CSV file
./sra-trawler.sh -c existing.csv -f reference.fasta

# Resume processing existing database
./sra-trawler.sh -d existing.db -f reference.fasta

# Resume and retry only download failures
./sra-trawler.sh -d existing.db -f reference.fasta -r

# Update existing database with new entries from CSV file
./sra-trawler.sh -c new_entries.csv -d existing.db -f reference.fasta

# Update existing database with new entries for a specific organism
./sra-trawler.sh -d existing.db -o "Aspergillus" -f reference.fasta
```

### Output
- `reads_mapping_to_reference/` - Compressed FASTQ files for entries with >1x coverage
- `alignments_mapping_to_reference/` - BAM alignment files for mapped reads
- `{organism}_sra_wgs.db` - SQLite database containing all SRA entries and status

In order to dump the database to a tsv file run:
` sqlite3 fungi_sra_wgs.db -header -separator $'\t' "SELECT * FROM sra_entries;" > fungi.tsv`

## 2. Plasmid Extractor (`plasmid_extractor.sh`)

The Plasmid Extractor script performs iterative assembly and refinement of plasmid sequences from sequencing reads.

### Purpose
- Takes sequencing reads and an initial reference plasmid
- Iteratively maps reads to the reference, assembles new contigs, and refines the plasmid sequence
- Uses both plasmid-specific and general assembly approaches
- Implements plateau detection to automatically terminate when no further improvement is possible
- Supports batch processing of multiple samples

### Limitation
- So far only processes paired illumina reads

### Usage Examples
```bash
# Process a single sample
./plasmid_extractor.sh --sample SRR123456 --initial-ref plasmid_reference.fasta

# Process all samples in a directory
./plasmid_extractor.sh --sample /path/to/reads/directory --batch --initial-ref plasmid_reference.fasta

# Adjust BLAST thresholds
./plasmid_extractor.sh --sample SRR123456 --initial-ref plasmid_reference.fasta --pid-threshold 80 --qcov-threshold 75

# Retry failed samples
./plasmid_extractor.sh --sample /path/to/reads/directory --batch --initial-ref plasmid_reference.fasta --retry
```

### Output
- `plasmids_extracted/plasmids/` - Final assembled plasmid sequences
- `plasmids_extracted/reads/` - Final round reads used for assembly
- `plasmids_extracted/working/` - Working directories for each sample
- Status tracking files for successful and failed samples

## 3. Plasmid Optimizer (`plasmid_optimizer.sh`)

The Plasmid Optimizer script post-processes extracted plasmids to improve their quality and consistency.

### Purpose
- Compares extracted plasmids against a reference plasmid
- Removes direct repeats that may have been introduced during assembly
- Adjusts start coordinates to match the reference
- Handles reverse complement orientation
- Optimizes plasmid sequences for downstream analysis

### Usage Examples
```bash
# Optimize plasmids in a directory
./plasmid_optimizer.sh --reference reference_plasmid.fasta --plasmids extracted_plasmids/

# Specify output directory and threads
./plasmid_optimizer.sh --reference reference_plasmid.fasta --plasmids extracted_plasmids/ --output optimized/
```

### Output
- Optimized plasmid sequences with terminal repeats removed
- Corrected orientations and start coordinates

# Installation

First install NCBI eDirect tools as described [here](https://www.nlm.nih.gov/dataguide/edirect/install.html).

For the rest of the tools needed the easiest way is to use the supplied `environment.yml` file and install it using conda or mamba.

`mamba env create -f environment.yml`

# Workflow

The typical workflow using these scripts is:

1. **SRA Trawler**: Search and download relevant SRA datasets
2. **Plasmid Extractor**: Extract and assemble plasmid sequences from the downloaded reads
3. **Plasmid Optimizer**: Optimize and standardize the extracted plasmids

This pipeline allows for systematic discovery and extraction of plasmid sequences from large-scale sequencing datasets.

