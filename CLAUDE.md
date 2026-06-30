# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Git Workflow Preferences

**IMPORTANT:** Do NOT automatically push to GitHub after making commits. Always:
1. Make code changes
2. Create local commits
3. **WAIT** for explicit user approval before pushing
4. Only run `git push` when the user explicitly requests it

This allows local testing before changes go to the remote repository.

## Overview

This repository contains a plasmid discovery pipeline that trawls through NCBI SRA databases to identify and extract plasmid sequences from sequencing reads. The pipeline consists of three main bash scripts that work sequentially or independently.

## Environment Setup

Install dependencies using conda/mamba:
```bash
mamba env create -f environment.yml
conda activate plasmidextractor
```

**Note:** NCBI eDirect tools must be installed separately following [NCBI documentation](https://www.nlm.nih.gov/dataguide/edirect/install.html).

## Core Scripts

### 1. sra-trawler.sh
Downloads and screens SRA entries for plasmid sequences. Manages a SQLite database to track processing status.

**Key functionality:**
- Downloads SRA entries via organism queries or CSV input
- Maps reads to reference genome using BWA-mem2 or minimap2
- Filters alignments to keep only primary alignments (excludes secondary/supplementary alignments with NNNNNN sequences)
- Retains only samples with >1x coverage
- Saves mapped reads and BAM alignments
- Tracks status in SQLite database with retry logic
- Optional push notifications when matches are found

**Common commands:**
```bash
# Create new database from organism query (default: fungi)
./sra-trawler.sh -f reference.fasta

# Resume processing existing database
./sra-trawler.sh -d existing.db -f reference.fasta

# Retry only failed downloads
./sra-trawler.sh -d existing.db -f reference.fasta -r

# Update database with new CSV entries
./sra-trawler.sh -c new_entries.csv -d existing.db -f reference.fasta

# Enable push notifications with custom message
./sra-trawler.sh -f reference.fasta -n "Plasmid" -m 5
```

**Database schema:**
```sql
CREATE TABLE sra_entries (
    run_accession TEXT PRIMARY KEY,
    sample_name TEXT,
    library_name TEXT,
    center_name TEXT,
    status TEXT,
    -- additional metadata fields
)
```

**Output structure:**
- `{organism}/{species}/reads/` - Compressed FASTQ files (>= min coverage only)
- `{organism}/{species}/alignments/` - BAM files for mapped reads
- `{organism}_sra_wgs.db` - SQLite database

### 2. plasmid_extractor.sh
Iteratively extracts and refines plasmid sequences through assembly-mapping cycles.

**Architecture - Iterative refinement loop:**
1. **Map reads to current reference** (BWA-mem2) → Extract mapped reads
2. **Paired-read integrity check** (seqkit) → Fix mismatched read counts
3. **Assembly attempt** (SPAdes):
   - First: plasmid mode (`--plasmid`)
   - Fallback: non-plasmid mode if plasmid mode fails
   - Subsample reads (50x, 30x, 20x coverage) if memory issues occur
4. **BLAST comparison** to initial reference → Filter contigs by PID/QCOV thresholds
5. **Select largest passing contig** as new reference for next round
6. **Convergence detection:**
   - Plateau: 3 consecutive rounds with same plasmid size (or <100bp growth)
   - No growth: 3 consecutive rounds with no read file size change
   - Safety: MAX_ROUNDS limit (default 100)

**Key environment variables:**
```bash
KMER_SET=21,33,55,77,99,127              # SPAdes k-mer sizes
MAX_ROUNDS=100                            # Maximum iteration rounds
BLAST_PID_THRESHOLD=70                    # Percent identity cutoff
BLAST_QCOV_THRESHOLD=70                   # Query coverage cutoff
RESULTS_DIR=plasmids_extracted           # Output directory
```

**Common commands:**
```bash
# Single sample
./plasmid_extractor.sh --sample SRR123456 --initial-ref plasmid_reference.fasta

# Batch processing directory
./plasmid_extractor.sh --sample /path/to/reads/ --batch --initial-ref plasmid_reference.fasta

# Adjust BLAST thresholds
./plasmid_extractor.sh --sample SRR123456 --initial-ref ref.fasta --pid-threshold 80 --qcov-threshold 75

# Retry failed samples
./plasmid_extractor.sh --sample /path/to/reads/ --batch --initial-ref ref.fasta --retry
```

**Input expectations:**
- Paired Illumina reads: `{sample}_1.fastq.gz` and `{sample}_2.fastq.gz`
- Limitation: Only processes paired Illumina reads currently

**Output structure:**
```
plasmids_extracted/
├── plasmids/           # Final assembled plasmid FASTA files
├── reads/              # Final round reads used for assembly
├── working/            # Per-sample working directories
│   └── {sample}/
│       ├── assembly/   # SPAdes output by round
│       ├── reads/      # Mapped reads by round
│       └── logs/       # Detailed logs including read_size_tracking.log
├── successful_samples.txt
└── failed_samples.txt
```

### 3. plasmid_optimizer.sh
Post-processes extracted plasmids to remove terminal repeats and standardize coordinates.

**Functionality:**
- BLASTs plasmids against reference
- Detects and removes direct terminal repeats
- Adjusts start coordinates to match reference
- Handles reverse complement orientation

**Common commands:**
```bash
# Optimize plasmids
./plasmid_optimizer.sh --reference reference.fasta --plasmids extracted_plasmids/

# Specify output and threads
./plasmid_optimizer.sh --reference ref.fasta --plasmids extracted_plasmids/ --output optimized/ --threads 8
```

### 4. local_mapper.sh
Maps local Illumina paired-end FASTQ files to a reference genome.

**Input format:** `xxx.R1.fastq.gz` and `xxx.R2.fastq.gz`

## Typical Workflow

```bash
# 1. Download and screen SRA datasets
./sra-trawler.sh -o "Aspergillus" -f reference_plasmid.fasta

# 2. Extract plasmids from downloaded reads
./plasmid_extractor.sh --sample aspergillus/aspergillus_fumigatus/reads/ --batch --initial-ref reference_plasmid.fasta

# 3. Optimize extracted plasmids
./plasmid_optimizer.sh --reference reference_plasmid.fasta --plasmids plasmids_extracted/plasmids/
```

## Architecture Details

### Convergence Strategy (plasmid_extractor.sh)
The extractor uses multiple convergence detection mechanisms:
- **Plateau detection**: Checks last 3 plasmid sizes; terminates if all equal or growth <100bp for 2 consecutive rounds
- **Read file tracking**: Monitors compressed FASTQ file sizes; terminates after 3 rounds of no change (unless good BLAST hit found)
- **Assembly fallback chain**: plasmid mode → non-plasmid mode → subsampled assembly (50x → 30x → 20x coverage)

### Database Management (sra-trawler.sh)
- Implements retry logic with exponential backoff for SQLite locking
- Supports resumption of interrupted processing
- Tracks per-sample status: pending, completed, failed, retry
- CSV format preserves sample names exactly as provided (e.g., "AR0037_AR0037")

### Tool Integration
- **Mapping**: BWA-mem2 (primary), minimap2 (fallback)
- **Assembly**: SPAdes with plasmid mode and k-mer optimization
- **QC**: fastp for read trimming, seqkit for paired-read validation
- **Analysis**: BLAST for contig filtering, samtools for BAM processing

## Debugging

### plasmid_extractor.sh logs:
- `logs/read_size_tracking.log` - Round-by-round read file size monitoring
- `logs/round{N}.spades_plasmid.log` - SPAdes plasmid mode output
- `logs/round{N}.spades_np.log` - SPAdes non-plasmid mode output
- `logs/round{N}.bwa_pe.log` - BWA paired-end mapping
- `logs/fastp.log` - Read trimming statistics

### sra-trawler.sh:
Enable debug mode: Set `DEBUG_MODE=true` at top of script

### Database queries:
```bash
# Export database to TSV
sqlite3 fungi_sra_wgs.db -header -separator $'\t' "SELECT * FROM sra_entries;" > fungi.tsv

# Check processing status
sqlite3 fungi_sra_wgs.db "SELECT status, COUNT(*) FROM sra_entries GROUP BY status;"
```

## Important Implementation Notes

### File Path Handling
Scripts use absolute paths internally. When providing input directories or files, relative paths are converted to absolute paths early in execution.

### Memory Management
plasmid_extractor.sh includes adaptive subsampling to handle large read sets that cause SPAdes to fail due to memory constraints.

### Read Integrity
The extractor validates paired-read counts after each mapping round and repairs mismatches using `seqkit pair` to prevent SPAdes errors.

### Hardcoded Paths
Recent commits mention fixing "hardcoded binary path" issues. When modifying scripts, ensure tool invocations use PATH-relative commands (e.g., `bwa-mem2` not `/usr/bin/bwa-mem2`) or check for tools in standard locations.
