#!/bin/bash

# Local Illumina Mapper Script
# This script maps local Illumina paired-end FASTQ files to a reference genome,
# and creates output folders containing the mapped reads and BAM files.
# 
# Input: Folder containing xxx.R1.fastq.gz and xxx.R2.fastq.gz files
# Output: Organized folders with mapped reads and BAM files

# Debug mode flag
DEBUG_MODE=false

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Global variable to track user interruption
INTERRUPTED=false

# Function to handle user interruption
handle_interrupt() {
    info_log "User interrupted (Ctrl+C). Stopping processing..."
    INTERRUPTED=true
    cleanup_on_exit
    exit 130  # Standard exit code for SIGINT
}

# Function to cleanup processes on exit
cleanup_on_exit() {
    info_log "Cleaning up processes..."
    
    # Kill any remaining mapping processes
    pkill -9 -f "bwa-mem2" 2>/dev/null
    pkill -9 -f "samtools" 2>/dev/null
}

# Function to print debug messages
debug_log() {
    :
}

# Function to print info messages
info_log() { echo "[INFO]  $*" >&2; }

# Function to print warning messages
warn_log() { :; }

# Function to print error messages
error_log() { echo "[ERROR] $*" >&2; }

# Function to print usage information
print_usage() {
    cat << EOF
Usage: $(basename "$0") [options]

Options:
    -h, --help          Show this help message and exit
    -i, --input DIR     Path to input directory containing FASTQ files (required)
                        Expects files named: xxx.R1.fastq.gz and xxx.R2.fastq.gz
    -r, --reference FILE Path to reference genome FASTA file (required)
                        Files will be mapped to this reference using BWA-MEM2
    -o, --output DIR    Output directory (default: ./mapped_output)
    -t, --threads INT   Number of threads for mapping (default: auto-detect)
    -c, --coverage FLOAT Minimum coverage threshold (default: 1.0)
                        Only samples with >= this coverage will be kept
    -D, --debug         Enable debug mode for verbose output

Examples:
    # Basic usage with default output directory and 1x coverage threshold
    $(basename "$0") -i /path/to/fastq/files -r reference.fasta
    
    # Specify custom output directory
    $(basename "$0") -i /path/to/fastq/files -r reference.fasta -o /path/to/output
    
    # Use specific number of threads and 5x coverage threshold
    $(basename "$0") -i /path/to/fastq/files -r reference.fasta -t 8 -c 5.0
    
    # Run with debug output enabled
    $(basename "$0") -i /path/to/fastq/files -r reference.fasta -D

Input Format:
    The input directory should contain paired-end FASTQ files with either naming convention:
    - sample1.R1.fastq.gz and sample1.R2.fastq.gz (dot separator)
    - sample1_R1.fastq.gz and sample1_R2.fastq.gz (underscore separator)
    - sample2.R1.fastq.gz and sample2.R2.fastq.gz
    - sample2_R1.fastq.gz and sample2_R2.fastq.gz
    etc.

Output Structure:
    {output_dir}/
    ├── reads/                    # Mapped reads (FASTQ) - only samples with sufficient coverage
    │   ├── sample1_R1.fastq.gz
    │   ├── sample1_R2.fastq.gz
    │   ├── sample2_R1.fastq.gz
    │   └── sample2_R2.fastq.gz
    ├── bam/                      # BAM alignment files - only samples with sufficient coverage
    │   ├── sample1.bam
    │   ├── sample1.bam.bai
    │   ├── sample2.bam
    │   └── sample2.bam.bai
    ├── logs/                     # Processing logs for all samples
    │   ├── sample1.log
    │   └── sample2.log
    └── coverage_results.txt      # Detailed coverage results for all samples
EOF
}

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to check required tools
check_required_tools() {
    local missing_tools=()
    
    if ! command_exists bwa-mem2; then
        missing_tools+=("bwa-mem2")
    fi
    
    if ! command_exists samtools; then
        missing_tools+=("samtools")
    fi
    
    if ! command_exists bc; then
        missing_tools+=("bc")
    fi
    
    if [[ ${#missing_tools[@]} -gt 0 ]]; then
        error_log "Missing required tools: ${missing_tools[*]}"
        error_log "Please install the missing tools and try again."
        exit 1
    fi
    
    info_log "All required tools are available"
}

# Function to find paired FASTQ files
find_paired_fastq_files() {
    local input_dir="$1"
    local pairs=()
    
    debug_log "Scanning directory: $input_dir"
    
    # Find all R1 files with different naming patterns
    local r1_patterns=("*.R1.fastq.gz" "*_R1.fastq.gz")
    
    for pattern in "${r1_patterns[@]}"; do
        debug_log "Searching for pattern: $pattern"
        
        while IFS= read -r -d '' r1_file; do
            local basename_r1=$(basename "$r1_file")
            local sample_name=""
            local r2_file=""
            
            # Extract sample name based on pattern
            if [[ "$pattern" == "*.R1.fastq.gz" ]]; then
                sample_name="${basename_r1%.R1.fastq.gz}"
                r2_file="${input_dir}/${sample_name}.R2.fastq.gz"
            elif [[ "$pattern" == "*_R1.fastq.gz" ]]; then
                sample_name="${basename_r1%_R1.fastq.gz}"
                r2_file="${input_dir}/${sample_name}_R2.fastq.gz"
            fi
            
            # Check if R2 file exists
            if [[ -f "$r2_file" ]]; then
                # Avoid duplicates
                local already_found=false
                for existing_pair in "${pairs[@]}"; do
                    local existing_sample=$(echo "$existing_pair" | cut -d: -f1)
                    if [[ "$existing_sample" == "$sample_name" ]]; then already_found=true; break; fi
                done
                if [[ "$already_found" == "false" ]]; then
                    pairs+=("$sample_name:$r1_file:$r2_file")
                fi
            else
                warn_log "No matching R2 file for: $r1_file (expected: $(basename "$r2_file"))"
            fi
        done < <(find "$input_dir" -maxdepth 1 -name "$pattern" -print0)
    done
    
    if [[ ${#pairs[@]} -eq 0 ]]; then
        error_log "No paired FASTQ files found in: $input_dir"
        exit 1
    fi
    
    # Return the pairs array by reference
    printf '%s\n' "${pairs[@]}"
}

# Function to calculate coverage (using the same approach as sra-trawler.sh)
calculate_coverage() {
    local bam_file="$1"; local reference="$2"; local threads="$3"; local redirect_output="$4"
    if [[ ! -f "$bam_file" ]]; then error_log "BAM file does not exist: $bam_file"; echo "0"; return; fi
    local read_count=$(samtools view -c "$bam_file" 2>/dev/null || echo "0")
    if [[ "$read_count" -eq 0 ]]; then error_log "BAM file has no reads"; echo "0"; return; fi
    [[ -f "${bam_file}.bai" ]] || samtools index "$bam_file" 1>/dev/null 2>/dev/null || { error_log "Failed to create BAM index"; echo "0"; return; }
    [[ -f "${reference}.fai" ]] || samtools faidx "$reference" 1>/dev/null 2>/dev/null || { error_log "Failed to create reference index"; echo "0"; return; }
    local ref_length=$(awk '{sum+=$2} END {print sum}' "${reference}.fai" 2>/dev/null || echo "0")
    if [[ "$ref_length" -eq 0 ]]; then error_log "Could not determine reference length"; echo "0"; return; fi
    local total_bases=$(samtools depth -a "$bam_file" 2>"$redirect_output" | awk '{sum += $3} END {print sum}')
    if [[ -z "$total_bases" ]] || [[ "$total_bases" == "0" ]]; then error_log "No depth data found or total bases is 0"; echo "0"; return; fi
    local coverage=$(echo "scale=2; $total_bases / $ref_length" | bc -l 2>/dev/null || echo "0")
    coverage="${coverage:-0.00}"; echo "$coverage"
}

# Function to process a single sample
process_sample() {
    local sample_name="$1"
    local r1_file="$2"
    local r2_file="$3"
    local reference="$4"
    local output_dir="$5"
    local threads="$6"
    local slot="$7"
    local min_coverage="${8:-1.0}"  # Default minimum coverage of 1x
    
    local temp_dir="${output_dir}/temp_${sample_name}"
    local log_file="${output_dir}/logs/${sample_name}.log"
    
    mkdir -p "$temp_dir"
    
    # Check if input files exist
    if [[ ! -f "$r1_file" ]]; then
        error_log "[Slot $slot] R1 file does not exist: $r1_file"
        return 1
    fi
    if [[ ! -f "$r2_file" ]]; then
        error_log "[Slot $slot] R2 file does not exist: $r2_file"
        return 1
    fi
    if [[ ! -f "$reference" ]]; then
        error_log "[Slot $slot] Reference file does not exist: $reference"
        return 1
    fi
    
    info_log "Processing: $sample_name"
    
    # Index reference if not already indexed
    if [[ ! -f "${reference}.amb" ]]; then
        if ! bwa-mem2 index "$reference" 1>/dev/null 2>/dev/null; then
            error_log "[Slot $slot] Failed to index reference genome"
            return 1
        fi
    fi

    # Map reads using BWA-MEM2 and create sorted BAM
    local redirect_output="/dev/null"
    [[ "$DEBUG_MODE" == "true" ]] && redirect_output="/dev/stderr"
    if ! bwa-mem2 mem -t "$threads" "$reference" "$r1_file" "$r2_file" 2>/dev/null | \
         samtools view -@ "$threads" -b -F 4 2>/dev/null | \
         samtools sort -@ "$threads" -o "${temp_dir}/${sample_name}_mapped.bam" 2>/dev/null; then
        error_log "[Slot $slot] Mapping failed for $sample_name"
        return 1
    fi

    # Create index for the BAM file
    if ! samtools index "${temp_dir}/${sample_name}_mapped.bam" 1>/dev/null 2>/dev/null; then
        error_log "[Slot $slot] Failed to create BAM index for $sample_name"
        return 1
    fi

    # Calculate coverage
    local coverage=$(calculate_coverage "${temp_dir}/${sample_name}_mapped.bam" "$reference" "$threads" "/dev/null")
    if [[ -z "$coverage" ]] || ! [[ "$coverage" =~ ^[0-9]+\.?[0-9]*$ ]]; then
        echo "FAILED"
        return 1
    fi

    # Compare to threshold
    local comparison_result=$(echo "$coverage >= $min_coverage" | bc -l 2>/dev/null || echo "0")
    
    # Handle success/failure outside the subshell
    if [[ "$comparison_result" == "1" ]]; then
        # Get mapping statistics directly from temp BAM
        local mapped_reads=$(samtools view -c "${temp_dir}/${sample_name}_mapped.bam" 2>/dev/null || echo "0")
        echo "SUCCESS:${coverage}:${mapped_reads}"
    else
        echo "FAILED:${coverage}:0"
    fi
    
    # Clean up temporary directory
    rm -rf "$temp_dir"
    
    return 0
}

# Function to process all samples
process_all_samples() {
    local reference="$1"
    local output_dir="$2"
    local threads="$3"
    local min_coverage="$4"
    shift 4  # Remove the first 4 parameters
    local -a pairs=("$@")  # Remaining parameters are the pairs
    
    local total_samples=${#pairs[@]}
    local passed=0
    local failed=0
    local coverage_failed=0
    
    info_log "Processing $total_samples samples sequentially"
    info_log "Minimum coverage threshold: ${min_coverage}x"
    
    # Create output directory (only for coverage files and summary)
    mkdir -p "${output_dir}"
    
    # Create coverage results file
    local coverage_file="${output_dir}/coverage_results.txt"
    echo -e "Sample\tStatus\tCoverage\tMapped_Reads" > "$coverage_file"
    
    # Process samples one by one
    for ((i=0; i<total_samples; i++)); do
        # Check for user interruption
        if [[ "$INTERRUPTED" == "true" ]]; then
            info_log "Processing interrupted by user. Stopping at sample $((i+1))/$total_samples"
            break
        fi
        
        local pair="${pairs[$i]}"
        local sample_name=$(echo "$pair" | cut -d: -f1)
        local r1_file=$(echo "$pair" | cut -d: -f2)
        local r2_file=$(echo "$pair" | cut -d: -f3)
        local slot=$((i + 1))
        
        info_log "Processing sample $slot/$total_samples: $sample_name"
        
        # Create temporary output file for this process
        local temp_output="${output_dir}/temp_${sample_name}_output.txt"
        
        # Run process_sample and capture output
        if process_sample "$sample_name" "$r1_file" "$r2_file" "$reference" "$output_dir" "$threads" "$slot" "$min_coverage" > "$temp_output"; then
            # Read the output to determine success/failure
            if [[ -f "$temp_output" ]]; then
                local result=$(cat "$temp_output")
                local status=$(echo "$result" | cut -d: -f1)
                local coverage=$(echo "$result" | cut -d: -f2)
                local mapped_reads=$(echo "$result" | cut -d: -f3)
                
                if [[ "$status" == "SUCCESS" ]]; then
                    ((passed++))
                    info_log "Sample $sample_name PASSED with ${coverage}x coverage ($passed/$total_samples)"
                    echo -e "${sample_name}\tPASSED\t${coverage}\t${mapped_reads}" >> "$coverage_file"
                else
                    ((coverage_failed++))
                    info_log "Sample $sample_name FAILED coverage check with ${coverage}x coverage ($coverage_failed/$total_samples)"
                    echo -e "${sample_name}\tFAILED_COVERAGE\t${coverage}\t0" >> "$coverage_file"
                fi
            else
                ((failed++))
                error_log "Sample $sample_name failed processing ($failed/$total_samples)"
                echo -e "${sample_name}\tFAILED_PROCESSING\t0\t0" >> "$coverage_file"
            fi
        else
            ((failed++))
            error_log "Sample $sample_name failed ($failed/$total_samples)"
            echo -e "${sample_name}\tFAILED_PROCESSING\t0\t0" >> "$coverage_file"
        fi
        
        # Clean up temporary output file
        rm -f "$temp_output"
        
        # Check for user interruption after each sample
        if [[ "$INTERRUPTED" == "true" ]]; then
            info_log "Processing interrupted by user. Stopping after sample $slot/$total_samples"
            break
        fi
    done
    
    # No separate summary file; coverage_results.txt is the single source of truth
    
    if [[ "$INTERRUPTED" == "true" ]]; then
        info_log "Processing interrupted by user. Completed: $passed passed, $coverage_failed failed coverage, $failed failed processing"
        return 130
    else
        info_log "Processing complete: $passed passed, $coverage_failed failed coverage, $failed failed processing"
        return $((failed + coverage_failed))
    fi
}

# Main function
main() {
    local input_dir=""
    local reference_file=""
    local output_dir="./mapped_output"
    local threads=$(nproc)
    local min_coverage=1.0
    
    # Parse command line arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                print_usage
                exit 0
                ;;
            -i|--input)
                input_dir="$2"
                shift 2
                ;;
            -r|--reference)
                reference_file="$2"
                shift 2
                ;;
            -o|--output)
                output_dir="$2"
                shift 2
                ;;
            -t|--threads)
                threads="$2"
                shift 2
                ;;
            -c|--coverage)
                min_coverage="$2"
                shift 2
                ;;
            -D|--debug)
                DEBUG_MODE=true
                shift
                ;;
            *)
                error_log "Unknown option: $1"
                print_usage
                exit 1
                ;;
        esac
    done
    
    # Validate required arguments
    if [[ -z "$input_dir" ]]; then
        error_log "Input directory is required"
        print_usage
        exit 1
    fi
    
    if [[ -z "$reference_file" ]]; then
        error_log "Reference file is required"
        print_usage
        exit 1
    fi
    
    if [[ ! -d "$input_dir" ]]; then
        error_log "Input directory does not exist: $input_dir"
        exit 1
    fi
    
    if [[ ! -f "$reference_file" ]]; then
        error_log "Reference file does not exist: $reference_file"
        exit 1
    fi
    
    # Check required tools
    check_required_tools
    
    # Setup trap to handle user interruption
    trap 'handle_interrupt' INT TERM
    
    info_log "Starting local Illumina mapper"
    info_log "Input directory: $input_dir"
    info_log "Reference file: $reference_file"
    info_log "Output directory: $output_dir"
    info_log "Threads: $threads"
    info_log "Minimum coverage threshold: ${min_coverage}x"
    
    # Find paired FASTQ files
    local pairs=()
    while IFS= read -r line; do
        if [[ -n "$line" ]]; then
            pairs+=("$line")
        fi
    done < <(find_paired_fastq_files "$input_dir")
    
    if [[ ${#pairs[@]} -eq 0 ]]; then
        error_log "No paired FASTQ files found"
        exit 1
    fi
    
    info_log "Found ${#pairs[@]} paired FASTQ files"
    
    # Process all samples
    process_all_samples "$reference_file" "$output_dir" "$threads" "$min_coverage" "${pairs[@]}"
    local exit_code=$?
    
    if [[ $exit_code -eq 0 ]]; then
        info_log "Completed"
    elif [[ $exit_code -eq 130 ]]; then
        info_log "Interrupted by user"
    else
        warn_log "Completed with errors"
    fi
    
    exit $exit_code
}

# Run the script
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
