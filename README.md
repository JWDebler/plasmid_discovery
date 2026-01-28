This is a selection of scripts I use to trawl through the NCBI SRA database using a reference sequence to find matches in sequencing reads that may not have made it into the respective assemblies.

# Scripts Overview

This repository contains four main scripts that work together to identify and extract plasmid sequences from sequencing data:

1. **`sra-trawler.sh`** - Downloads and screens SRA entries for potential plasmid (or other) sequences
2. **`plasmid_extractor.sh`** - Iteratively extracts and refines plasmid sequences from sequencing reads
3. **`plasmid_optimizer.sh`** - Optimizes extracted plasmids by removing repeats and adjusting coordinates
4. **`local_mapper.sh`** - Maps local Illumina paired-end reads to a reference genome

## 1. SRA Trawler (`sra-trawler.sh`)

The SRA Trawler script is designed to systematically search through the NCBI SRA database to find sequencing datasets that contain potential plasmid sequences matching a supplied reference.

### Purpose
- Downloads all SRA entries for a given organism queries or CSV input
- Maps downloaded sequencing reads to a reference genome
- Filters alignments to keep only primary alignments (excludes secondary/supplementary alignments)
- Identifies entries with sufficient coverage (>1x) as potential plasmid sources
- Saves only reads that map to the reference for downstream analysis
- Maintains a SQLite database to track processing status
- Optional push notifications when matches are found (requires [pushbullet-bash](https://github.com/Red5d/pushbullet-bash))


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

# Enable push notifications when matches are found
./sra-trawler.sh -f reference.fasta -n "Plasmid"

# Adjust minimum coverage threshold
./sra-trawler.sh -f reference.fasta -m 5

# Combine multiple options
./sra-trawler.sh -o "Aspergillus" -f reference.fasta -m 5 -n "Plasmid" -x 10
```

### Key Options
- `-f, --reference FILE` - Reference genome FASTA file (required)
- `-o, --organism STR` - Organism to search for in SRA (default: fungi)
- `-c, --csv FILE` - Input CSV file with SRA accessions
- `-d, --db FILE` - Path to SQLite database file
- `-m, --min-coverage NUM` - Minimum coverage threshold for saving reads (default: 1)
- `-n, --notify STR` - Enable push notifications with custom message string (e.g., "Plasmid" will notify "Plasmid detected! ...")
- `-x, --connections INT` - Number of concurrent connections to SRA (default: 5)
- `-r, --retry` - Retry failed downloads only (selective retry)
- `-R, --retry-all` - Retry ALL failed entries (comprehensive retry)
- `-D, --debug` - Enable debug mode for verbose output

### Push Notifications
When the `-n/--notify` option is used, the script will send push notifications via Pushbullet when samples meeting the coverage threshold are detected. For example:

```bash
./sra-trawler.sh -f plasmid.fasta -n "Plasmid"
```

Will send notifications like: `Plasmid detected! SRR123456 (Aspergillus fumigatus) has 5.2x coverage`

**Requirements:**
- Install [pushbullet-bash](https://github.com/Red5d/pushbullet-bash) (`pb` command)
- Configure Pushbullet API token as described in the tool's documentation
- The script will check if `pb` is installed and exit with an error if not found

### Output
The script creates an organized directory structure based on organism and species:
```
{organism}/                          # Top-level folder (e.g., "fungi", "aspergillus")
└── {species}/                       # Species subfolder (e.g., "aspergillus_fumigatus")
    ├── reads/                       # Compressed FASTQ files (>= min coverage only)
    │   ├── SRR123456_1.fastq.gz
    │   └── SRR123456_2.fastq.gz
    └── alignments/                  # BAM alignment files for mapped reads
        └── SRR123456_mapped.bam

{organism}_sra_wgs.db                # SQLite database with all entries and status
```

**Database export:**
To dump the database to a TSV file:
```bash
sqlite3 fungi_sra_wgs.db -header -separator $'\t' "SELECT * FROM sra_entries;" > fungi.tsv
```

**Check processing status:**
```bash
sqlite3 fungi_sra_wgs.db "SELECT status, COUNT(*) FROM sra_entries GROUP BY status;"
```

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

## 4. Local Mapper (`local_mapper.sh`)

The Local Mapper script maps local Illumina paired-end FASTQ files to a reference genome.

### Purpose
- Maps local sequencing reads to a reference genome
- Handles paired-end Illumina data
- Useful for analyzing your own sequencing data against a reference

### Input Format
- Expects paired-end reads with naming pattern: `xxx.R1.fastq.gz` and `xxx.R2.fastq.gz`

# Installation

First install NCBI eDirect tools as described [here](https://www.nlm.nih.gov/dataguide/edirect/install.html).

For the rest of the tools needed the easiest way is to use the supplied `environment.yml` file and install it using conda or mamba.

`mamba env create -f environment.yml`

# Typical Workflow

The typical workflow using these scripts is:

## 1. Discovery Phase (SRA Trawler)
Search and download relevant SRA datasets that contain your plasmid of interest:

```bash
# Start with a reference plasmid sequence
./sra-trawler.sh -o "Aspergillus" -f reference_plasmid.fasta -m 5 -n "Plasmid"
```

This will:
- Query SRA for all Aspergillus WGS datasets
- Map reads to your reference plasmid
- Keep only samples with ≥5x coverage
- Send push notifications when matches are found
- Save reads and alignments organized by species

## 2. Extraction Phase (Plasmid Extractor)
Extract and assemble complete plasmid sequences from the downloaded reads:

```bash
# Process all downloaded samples in batch mode
./plasmid_extractor.sh --sample fungi/aspergillus_fumigatus/reads/ --batch --initial-ref reference_plasmid.fasta
```

This will iteratively refine each plasmid through multiple assembly rounds until convergence.

## 3. Optimization Phase (Plasmid Optimizer)
Post-process the extracted plasmids to remove artifacts and standardize coordinates:

```bash
# Optimize all extracted plasmids
./plasmid_optimizer.sh --reference reference_plasmid.fasta --plasmids plasmids_extracted/plasmids/
```

This pipeline allows for systematic discovery and extraction of plasmid sequences from large-scale sequencing datasets.

