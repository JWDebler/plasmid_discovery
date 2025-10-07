#!/usr/bin/env bash
# Temporarily disable strict error handling to allow graceful handling of errors
# set -euo pipefail



SPADES_BIN=${SPADES_BIN:-/opt/spades/SPAdes-4.2.0-Linux/bin/spades.py}
KMER_SET=${KMER_SET:-21,33,55,77,99,127}
MAX_ROUNDS=${MAX_ROUNDS:-100}

# BLAST thresholds (override via environment variables)
BLAST_PID_THRESHOLD=${BLAST_PID_THRESHOLD:-70}
BLAST_QCOV_THRESHOLD=${BLAST_QCOV_THRESHOLD:-70}

usage() {
        echo "Usage: $0 [--sample <sample_id_or_dir>] [--initial-ref <ref.fa>] [--max-rounds <N>] [--batch]" >&2
        echo "  --sample, -s        Sample ID (prefix before _1.fastq.gz) OR directory containing *_1.fastq.gz and *_2.fastq.gz" >&2
        echo "  --initial-ref, -r   Initial reference FASTA" >&2
        echo "  --max-rounds, -m    Maximum iteration rounds (default: 100 or env MAX_ROUNDS)" >&2
        echo "  --batch, -b         Process all sample pairs in directory (use with --sample <dir>)" >&2
        echo "  --output-dir, -o    Output base directory (default: plasmids_extracted or env RESULTS_DIR)" >&2
        echo "  --pid-threshold     BLAST percent identity cutoff (default: ${BLAST_PID_THRESHOLD})" >&2
        echo "  --qcov-threshold    BLAST query coverage cutoff (default: ${BLAST_QCOV_THRESHOLD})" >&2
        echo "  --retry             Reprocess samples listed as failed (ignores failure cache)" >&2
        echo "  --help, -h          Show this help" >&2
}

derive_sample_from_dir() {
        local dir="$1"
        local first_base=""
        shopt -s nullglob
        for f in "$dir"/*_1.fastq.gz; do
                local name
                name=$(basename "$f")
                local base=${name%_1.fastq.gz}
                if [ -f "$dir/${base}_2.fastq.gz" ]; then
                        first_base="$base"
                        break
                fi
        done
        shopt -u nullglob
        if [ -z "$first_base" ]; then
                echo "No paired *_1.fastq.gz and *_2.fastq.gz found in: $dir" >&2
                exit 1
        fi
        echo "$first_base"
}

get_all_samples_from_dir() {
        local dir="$1"
        local samples=()
        shopt -s nullglob
        echo "Scanning directory: $dir" >&2
        for f in "$dir"/*_1.fastq.gz; do
                local name
                name=$(basename "$f")
                local base=${name%_1.fastq.gz}
                echo "Found *_1.fastq.gz: $name -> base: $base" >&2
                if [ -f "$dir/${base}_2.fastq.gz" ]; then
                        echo "  Paired with: ${base}_2.fastq.gz ✓" >&2
                        samples+=("$base")
                else
                        echo "  No matching *_2.fastq.gz found ✗" >&2
                fi
        done
        shopt -u nullglob
        if [ ${#samples[@]} -eq 0 ]; then
                echo "No paired *_1.fastq.gz and *_2.fastq.gz found in: $dir" >&2
                exit 1
        fi
        echo "Total samples found: ${#samples[@]}" >&2
        printf '%s\n' "${samples[@]}"
}

# Parse CLI args
SAMPLE_ARG=""
INITIAL_REF_ARG=""
MAX_ROUNDS_ARG=""
BATCH_MODE=false
# Retry mode: ignore previously failed samples and try again
RETRY_MODE=false
OUTPUT_DIR_ARG=""
while [[ $# -gt 0 ]]; do
        case "$1" in
                --sample|-s)
                        SAMPLE_ARG="${2:-}"
                        shift 2
                        ;;
                --initial-ref|-r)
                        INITIAL_REF_ARG="${2:-}"
                        shift 2
                        ;;
                --max-rounds|-m)
                        MAX_ROUNDS_ARG="${2:-}"
                        shift 2
                        ;;
                --batch|-b)
                        BATCH_MODE=true
                        shift
                        ;;
                --output-dir|-o)
                        OUTPUT_DIR_ARG="${2:-}"
                        shift 2
                        ;;
                --pid-threshold)
                        BLAST_PID_THRESHOLD="${2:-}"
                        shift 2
                        ;;
                --qcov-threshold)
                        BLAST_QCOV_THRESHOLD="${2:-}"
                        shift 2
                        ;;
                --retry)
                        RETRY_MODE=true
                        shift
                        ;;
                --help|-h)
                        usage
                        exit 0
                        ;;
                *)
                        echo "Unknown argument: $1" >&2
                        usage
                        exit 1
                        ;;
        esac
done

# Resolve config from args/defaults
RAW_DIR="."
BATCH_SAMPLES=()
# Base results directory can be provided via --output-dir, else fallback to env RESULTS_DIR, else default
RESULTS_DIR_DEFAULT=${RESULTS_DIR:-"plasmids_extracted"}
if [ -n "${OUTPUT_DIR_ARG}" ]; then
        RESULTS_DIR="${OUTPUT_DIR_ARG}"
else
        RESULTS_DIR="${RESULTS_DIR_DEFAULT}"
fi
WORKING_DIR="${RESULTS_DIR}/working"  # New working directory for individual samples
PLASMIDS_DIR="${RESULTS_DIR}/plasmids"  # Directory for final plasmids
READS_DIR="${RESULTS_DIR}/reads"  # Directory for final reads
mkdir -p "${RESULTS_DIR}" "${WORKING_DIR}" "${PLASMIDS_DIR}" "${READS_DIR}"

# Initialize tracking files
FAILED_SAMPLES_FILE="${RESULTS_DIR}/failed_samples.txt"
SUCCESS_SAMPLES_FILE="${RESULTS_DIR}/success_samples.txt"
PROCESSING_STATUS_FILE="${RESULTS_DIR}/processing_status.txt"

# Check for existing status files and report on already processed samples
echo "Checking for existing processing status..." >&2
if [ -f "${SUCCESS_SAMPLES_FILE}" ] && [ -s "${SUCCESS_SAMPLES_FILE}" ]; then
        success_count=$(wc -l < "${SUCCESS_SAMPLES_FILE}")
        echo "Found ${success_count} previously successfully processed samples" >&2
fi
if [ -f "${FAILED_SAMPLES_FILE}" ] && [ -s "${FAILED_SAMPLES_FILE}" ]; then
        failed_count=$(wc -l < "${FAILED_SAMPLES_FILE}")
        echo "Found ${failed_count} previously failed samples" >&2
fi

# Create tracking files if they don't exist (don't clear existing ones)
touch "${FAILED_SAMPLES_FILE}" "${SUCCESS_SAMPLES_FILE}" "${PROCESSING_STATUS_FILE}"

if [ -n "$SAMPLE_ARG" ]; then
        if [ -d "$SAMPLE_ARG" ]; then
            RAW_DIR="$(cd "$SAMPLE_ARG" && pwd)"
            echo "Input is a directory: ${RAW_DIR}" >&2
            if [ "$BATCH_MODE" = true ]; then
                echo "Batch mode: processing all samples in directory: ${RAW_DIR}" >&2
                mapfile -t BATCH_SAMPLES < <(get_all_samples_from_dir "$RAW_DIR")
                echo "Found ${#BATCH_SAMPLES[@]} samples: ${BATCH_SAMPLES[*]}" >&2
            else
                echo "Single sample mode: detecting first sample from directory: ${RAW_DIR}" >&2
                SAMPLE="$(derive_sample_from_dir "$RAW_DIR")"
                echo "Detected sample '${SAMPLE}' from directory: ${RAW_DIR}" >&2
                BATCH_SAMPLES=("$SAMPLE")
            fi
        else
            # Check if the argument contains a path (e.g., "reads/SRR32106019")
            if [[ "$SAMPLE_ARG" == */* ]]; then
                # Extract directory and sample ID from path
                sample_dir="$(dirname "$SAMPLE_ARG")"
                sample_id="$(basename "$SAMPLE_ARG")"
                
                # Check if the directory exists
                if [ -d "$sample_dir" ]; then
                    RAW_DIR="$(cd "$sample_dir" && pwd)"
                    SAMPLE="$sample_id"
                    echo "Input is a single sample ID with directory path: ${SAMPLE_ARG}" >&2
                    echo "Extracted directory: ${RAW_DIR}" >&2
                    echo "Extracted sample ID: ${SAMPLE}" >&2
                    BATCH_SAMPLES=("$SAMPLE")
                else
                    echo "ERROR: Directory '${sample_dir}' does not exist for sample path '${SAMPLE_ARG}'" >&2
                    exit 1
                fi
            else
                echo "Input is a single sample ID: ${SAMPLE_ARG}" >&2
                SAMPLE="$SAMPLE_ARG"
                BATCH_SAMPLES=("$SAMPLE")
            fi
        fi
else
        echo "ERROR: No sample specified. Please provide a sample using --sample or -s option." >&2
        echo "Use --help for usage information." >&2
        exit 1
fi

echo "Final BATCH_SAMPLES array: [${BATCH_SAMPLES[*]}]" >&2
echo "BATCH_MODE: ${BATCH_MODE}" >&2
echo "RETRY_MODE: ${RETRY_MODE}" >&2
echo "BLAST thresholds -> PID: ${BLAST_PID_THRESHOLD}, QCOV: ${BLAST_QCOV_THRESHOLD}" >&2

if [ -n "$INITIAL_REF_ARG" ]; then
        INITIAL_REF="$INITIAL_REF_ARG"
else
        INITIAL_REF="$INITIAL_REF_DEFAULT"
fi

if [ -n "$MAX_ROUNDS_ARG" ]; then
        MAX_ROUNDS="$MAX_ROUNDS_ARG"
fi

# Function to safely get seqkit stats (handles header line)
get_seqkit_count() {
        local file="$1"
        local count=$(seqkit stats -T "$file" 2>/dev/null | grep -v "^file" | tail -n1 | cut -f4)
        # Validate that it's a number
        if [[ "$count" =~ ^[0-9]+$ ]]; then
                echo "$count"
        else
                echo "0"
        fi
}

# Function to safely get seqkit length (handles header line)
get_seqkit_length() {
        local file="$1"
        local length=$(seqkit stats -T "$file" 2>/dev/null | grep -v "^file" | tail -n1 | cut -f2)
        # Validate that it's a number
        if [[ "$length" =~ ^[0-9]+$ ]]; then
                echo "$length"
        else
                echo "0"
        fi
}

# Function to subsample reads when assembly fails
subsample_reads_for_assembly() {
        local reads_dir="$1"
        local round="$2"
        local target_coverage="${3:-50}"  # Default 50x coverage
        local sample_name="$4"
        
        echo "Round ${round}: Assembly failed. Attempting read subsampling with ${target_coverage}x target coverage..." >&2
        
        # Calculate current read counts (safely)
        local r1_count=$(get_seqkit_count "${reads_dir}/round${round}.1.fastq.gz")
        local r2_count=$(get_seqkit_count "${reads_dir}/round${round}.2.fastq.gz")
        local s_count=$(get_seqkit_count "${reads_dir}/round${round}.s.fastq.gz")
        
        echo "Round ${round}: Current reads - R1: ${r1_count}, R2: ${r2_count}, S: ${s_count}" >&2
        
        # Estimate current coverage (rough calculation)
        local current_ref_length=0
        if [ -f "${current_ref}" ]; then
                current_ref_length=$(get_seqkit_length "${current_ref}")
        fi
        
        if [ "$current_ref_length" -gt 0 ]; then
                local current_coverage=$(echo "scale=1; (${r1_count} + ${r2_count} + ${s_count}) * 150 / ${current_ref_length}" | bc -l)
                echo "Round ${round}: Estimated current coverage: ${current_coverage}x" >&2
                
                # Calculate subsampling ratio
                local subsample_ratio=$(echo "scale=3; ${target_coverage} / ${current_coverage}" | bc -l)
                
                # Ensure ratio is between 0.1 and 1.0
                if (( $(echo "$subsample_ratio < 0.1" | bc -l) )); then
                        subsample_ratio=0.1
                elif (( $(echo "$subsample_ratio > 1.0" | bc -l) )); then
                        subsample_ratio=1.0
                fi
                
                echo "Round ${round}: Subsampling ratio: ${subsample_ratio}" >&2
                
                # Create subsampled read files
                local subsampled_dir="${reads_dir}/subsampled_round${round}_${target_coverage}x"
                mkdir -p "${subsampled_dir}"
                
                # Subsample paired reads
                seqkit sample -p "${subsample_ratio}" "${reads_dir}/round${round}.1.fastq.gz" -o "${subsampled_dir}/round${round}.1.fastq.gz" 2>/dev/null
                seqkit sample -p "${subsample_ratio}" "${reads_dir}/round${round}.2.fastq.gz" -o "${subsampled_dir}/round${round}.2.fastq.gz" 2>/dev/null
                
                # Subsample singleton reads
                seqkit sample -p "${subsample_ratio}" "${reads_dir}/round${round}.s.fastq.gz" -o "${subsampled_dir}/round${round}.s.fastq.gz" 2>/dev/null
                
                # Check subsampled read counts (safely)
                local sub_r1_count=$(get_seqkit_count "${subsampled_dir}/round${round}.1.fastq.gz")
                local sub_r2_count=$(get_seqkit_count "${subsampled_dir}/round${round}.2.fastq.gz")
                local sub_s_count=$(get_seqkit_count "${subsampled_dir}/round${round}.s.fastq.gz")
                
                echo "Round ${round}: Subsampled reads - R1: ${sub_r1_count}, R2: ${sub_r2_count}, S: ${sub_s_count}" >&2
                
                # Return the subsampled directory path
                echo "${subsampled_dir}"
        else
                echo "Round ${round}: Cannot estimate coverage, using original reads" >&2
                echo "${reads_dir}"
        fi
}

# Function to check required tools
check_required_tools() {
        local missing_tools=()
        
        # Check for required tools
        for tool in bwa-mem2 samtools fastp seqkit pigz; do
                if ! command -v "$tool" >/dev/null 2>&1; then
                        missing_tools+=("$tool")
                fi
        done
        
        # Check for SPAdes
        if [ ! -f "${SPADES_BIN}" ] && ! command -v spades.py >/dev/null 2>&1; then
                missing_tools+=("spades.py")
        fi
        
        if [ ${#missing_tools[@]} -gt 0 ]; then
                echo "ERROR: Missing required tools: ${missing_tools[*]}" >&2
                echo "Please install the missing tools and try again." >&2
                return 1
        fi
        
        echo "All required tools are available ✓" >&2
        return 0
}

# Function to check if sample was already processed
is_sample_processed() {
        local sample="$1"
        if grep -q "^${sample}$" "${SUCCESS_SAMPLES_FILE}" 2>/dev/null; then
                return 0  # Sample was successfully processed
        fi
        # If retry mode is enabled, ignore failed cache
        if [ "${RETRY_MODE}" != true ]; then
                if grep -q "^${sample}$" "${FAILED_SAMPLES_FILE}" 2>/dev/null; then
                        return 0  # Sample failed previously
                fi
        fi
        return 1  # Sample not processed yet
}

# Function to mark sample as processed
mark_sample_processed() {
        local sample="$1"
        local status="$2"  # "success" or "failed"
        local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
        
        if [ "$status" = "success" ]; then
                echo "${sample}" >> "${SUCCESS_SAMPLES_FILE}"
                printf "%s\tSUCCESS\t%s\n" "${sample}" "${timestamp}" >> "${PROCESSING_STATUS_FILE}"
        else
                echo "${sample}" >> "${FAILED_SAMPLES_FILE}"
                printf "%s\tFAILED\t%s\n" "${sample}" "${timestamp}" >> "${PROCESSING_STATUS_FILE}"
        fi
}

# Function to detect plateau in plasmid size growth
detect_plateau() {
        local size_history=("$@")
        local history_size=${#size_history[@]}
        
        # Need at least 3 data points to detect plateau
        if [ $history_size -lt 3 ]; then
                return 1  # No plateau detected
        fi
        
        # Get the last 3 sizes
        local last_3_sizes=()
        for ((i=history_size-3; i<history_size; i++)); do
                last_3_sizes+=("${size_history[$i]}")
        done
        
        # Check if all last 3 sizes are the same (exact plateau)
        if [ "${last_3_sizes[0]}" -eq "${last_3_sizes[1]}" ] && [ "${last_3_sizes[1]}" -eq "${last_3_sizes[2]}" ]; then
                return 0  # Plateau detected
        fi
        
        # Check if growth has slowed significantly (less than 100 bp growth for 2 consecutive rounds)
        local growth_1=$((last_3_sizes[1] - last_3_sizes[0]))
        local growth_2=$((last_3_sizes[2] - last_3_sizes[1]))
        
        # If growth is very small (< 100 bp) for 2 consecutive rounds, consider it a plateau
        if [ ${growth_1} -lt 100 ] && [ ${growth_2} -lt 100 ]; then
                return 0  # Plateau detected
        fi
        
        return 1  # No plateau detected
}

process_sample() {
        local sample="$1"
        
        echo "=== Sample: ${sample}" >&2
        echo "======================================================" >&2
        
        # Check if sample was already processed
        if is_sample_processed "$sample"; then
                echo "Sample ${sample} was already processed. Skipping." >&2
                return 0
        fi
        
        echo "Sample ${sample} not previously processed. Starting processing..." >&2
        
        local sample_dir="${WORKING_DIR}/${sample}"  # Use working directory
        local log_dir="${sample_dir}/logs"
        local asm_dir="${sample_dir}/assembly"
        local reads_dir="${sample_dir}/reads"
        mkdir -p "${log_dir}" "${asm_dir}" "${reads_dir}"
        
        local r1_raw="${RAW_DIR}/${sample}_1.fastq.gz"
        local r2_raw="${RAW_DIR}/${sample}_2.fastq.gz"
        
        # Prepare BLAST DB for the initial reference (once per sample)
        if [ -f "${INITIAL_REF}" ]; then
                makeblastdb -in "${INITIAL_REF}" -dbtype nucl >/dev/null 2>&1 || true
                # Get the base name of the reference file for BLAST database (keep .fasta part)
                BLAST_DB_BASE="${INITIAL_REF}"
        else
                echo "Initial reference not found: ${INITIAL_REF}" >&2
                mark_sample_processed "$sample" "failed"
                return 0  # Return 0 to continue with next sample instead of terminating
        fi
        
        # Check if trimmed reads already exist
        local r1_trim="${reads_dir}/${sample}_1.trimmed.fastq.gz"
        local r2_trim="${reads_dir}/${sample}_2.trimmed.fastq.gz"
        local single_trim="${reads_dir}/${sample}_single.trimmed.fastq.gz"
        
        if [ -f "${r1_trim}" ] && [ -f "${r2_trim}" ] && [ -f "${single_trim}" ]; then
                echo "Trimmed reads already exist for sample ${sample}. Skipping trimming step."
                echo "Using existing trimmed reads:"
                echo "  - ${r1_trim}"
                echo "  - ${r2_trim}"
                echo "  - ${single_trim}"
        else
                echo "Trimming reads for sample ${sample}..."
                {
                fastp \
                        --in1 "${r1_raw}" \
                        --in2 "${r2_raw}" \
                        --out1 "${r1_trim}" \
                        --out2 "${r2_trim}" \
                        --unpaired1 "${reads_dir}/${sample}_single1.trimmed.fastq.gz" \
                        --unpaired2 "${reads_dir}/${sample}_single2.trimmed.fastq.gz" \
                        --trim_poly_g \
                        --cut_right \
                        --thread 16 \
                        --disable_quality_filtering
                } > "${log_dir}/fastp.log" 2>&1 &
                local fp_pid=$!
                echo "fastp launched (PID ${fp_pid}). Waiting..."
                wait ${fp_pid}
                echo "Trimming complete. Logs: ${log_dir}/fastp.log"
                
                # Combine single reads
                cat "${reads_dir}/${sample}_single1.trimmed.fastq.gz" "${reads_dir}/${sample}_single2.trimmed.fastq.gz" > "${single_trim}" || true
                rm -f "${reads_dir}/${sample}_single1.trimmed.fastq.gz" "${reads_dir}/${sample}_single2.trimmed.fastq.gz" || true
        fi
        
        local current_ref="${INITIAL_REF}"
        local prev_max_len=0
        local round=1
        local largest_fa=""
        local found_good_blast_hit=false  # Track if we found any good BLAST hit
        local found_good_np_blast_hit=false  # Track if a good BLAST hit came from non-plasmid SPAdes
        local plateau_detected=false  # Track if we exited due to plateau detection
        
        # NEW: Read file size tracking variables
        local prev_read_size=0
        local no_growth_cycles=0
        local max_no_growth_cycles=3
        local size_tracking_log="${log_dir}/read_size_tracking.log"
        echo "Round\tR1_Size\tR2_Size\tS_Size\tTotal_Size\tGrowth\tNo_Growth_Cycles" > "${size_tracking_log}"
        
        # NEW: Plasmid size history for plateau detection
        local plasmid_size_history=()
        
        while [ ${round} -le ${MAX_ROUNDS} ]; do
                local outdir="${asm_dir}/round${round}"
                mkdir -p "${outdir}"
                
                # Mapping original reads to current reference and extracting mapped reads
                echo "Round ${round}: Indexing reference and mapping paired reads to ${current_ref}..."
                bwa-mem2 index "${current_ref}" >/dev/null 2>&1 || true
                { bwa-mem2 mem -t $(nproc) "${current_ref}" "${r1_trim}" "${r2_trim}" \
                        | samtools view -@ $(nproc) -b -F 4 \
                        | samtools fastq -@ $(nproc) \
                                -1 "${reads_dir}/round${round}.1.fastq.gz" \
                                -2 "${reads_dir}/round${round}.2.fastq.gz" \
                                -s "${reads_dir}/round${round}.s1.fastq.gz" \
                                -0 /dev/null -n; } 2> "${log_dir}/round${round}.bwa_pe.log"
                echo "Round ${round}:    mapping paired reads complete."
                
                # NEW: Sanity check and fix paired read files using seqkit
                echo "Round ${round}:    checking paired read file integrity..."
                local r1_file="${reads_dir}/round${round}.1.fastq.gz"
                local r2_file="${reads_dir}/round${round}.2.fastq.gz"
                
                if [ -f "${r1_file}" ] && [ -f "${r2_file}" ]; then
                    # Check if files have unequal number of reads
                    local r1_count=$(seqkit stats -T "${r1_file}" | tail -n1 | cut -f4)
                    local r2_count=$(seqkit stats -T "${r2_file}" | tail -n1 | cut -f4)
                    
                    echo "Round ${round}:    R1 reads: ${r1_count}, R2 reads: ${r2_count}" >&2
                    
                    if [ "${r1_count}" != "${r2_count}" ]; then
                        echo "Round ${round}:    WARNING - Unequal read counts detected! Fixing with seqkit..." >&2
                        
                        # Create backup of original files
                        cp "${r1_file}" "${r1_file}.backup"
                        cp "${r2_file}" "${r2_file}.backup"
                        
                        # Use seqkit to fix paired reads (this will ensure equal counts)
                        seqkit pair -1 "${r1_file}.backup" -2 "${r2_file}.backup" -o "${reads_dir}/round${round}.fixed" 2>"${log_dir}/round${round}.seqkit_fix.log"

                        # Prefer explicit fixed outputs; otherwise handle seqkit-created *.paired.backup cases
                        if [ -f "${reads_dir}/round${round}.fixed.1.fastq.gz" ] && [ -f "${reads_dir}/round${round}.fixed.2.fastq.gz" ]; then
                            mv -f "${reads_dir}/round${round}.fixed.1.fastq.gz" "${r1_file}"
                            mv -f "${reads_dir}/round${round}.fixed.2.fastq.gz" "${r2_file}"
                        elif [ -f "${r1_file}.paired.backup" ] && [ -f "${r2_file}.paired.backup" ]; then
                            # Some environments/tools generate *.paired.backup as the corrected files
                            mv -f "${r1_file}.paired.backup" "${r1_file}"
                            mv -f "${r2_file}.paired.backup" "${r2_file}"
                        else
                            echo "Round ${round}:    ERROR - Expected fixed paired files not found. Restoring originals." >&2
                            mv -f "${r1_file}.backup" "${r1_file}"
                            mv -f "${r2_file}.backup" "${r2_file}"
                        fi

                        # Check the fixed counts
                        local fixed_r1_count=$(seqkit stats -T "${r1_file}" | tail -n1 | cut -f4)
                        local fixed_r2_count=$(seqkit stats -T "${r2_file}" | tail -n1 | cut -f4)
                        echo "Round ${round}:    Fixed! R1 reads: ${fixed_r1_count}, R2 reads: ${fixed_r2_count}" >&2
                        
                        # Clean up backup and intermediate files if present
                        rm -f "${r1_file}.backup" "${r2_file}.backup" \
                              "${r1_file}.paired.backup" "${r2_file}.paired.backup" \
                              "${reads_dir}/round${round}.fixed.1.fastq.gz" "${reads_dir}/round${round}.fixed.2.fastq.gz"
                    else
                        echo "Round ${round}:    Paired read files have equal counts ✓" >&2
                    fi
                else
                    echo "Round ${round}:    WARNING - One or both paired read files missing!" >&2
                fi
                
                echo "Round ${round}: Mapping singleton reads to ${current_ref}..."
                { bwa-mem2 mem -t $(nproc) "${current_ref}" "${single_trim}" \
                        | samtools view -@ $(nproc) -b -F 4 \
                        | samtools fastq -@ $(nproc) | pigz > "${reads_dir}/round${round}.s2.fastq.gz"; } 2> "${log_dir}/round${round}.bwa_se.log"
                echo "Round ${round}:    mapping singleton reads complete."
                
                cat "${reads_dir}/round${round}.s1.fastq.gz" "${reads_dir}/round${round}.s2.fastq.gz" > "${reads_dir}/round${round}.s.fastq.gz"
                rm -f "${reads_dir}/round${round}.s1.fastq.gz" "${reads_dir}/round${round}.s2.fastq.gz"
                
                # NEW: Calculate current read file sizes after mapping
                local r1_size=0
                local r2_size=0
                local s_size=0
                local total_read_size=0
                
                if [ -f "${reads_dir}/round${round}.1.fastq.gz" ]; then
                        r1_size=$(stat -c%s "${reads_dir}/round${round}.1.fastq.gz" 2>/dev/null || echo 0)
                fi
                if [ -f "${reads_dir}/round${round}.2.fastq.gz" ]; then
                        r2_size=$(stat -c%s "${reads_dir}/round${round}.2.fastq.gz" 2>/dev/null || echo 0)
                fi
                if [ -f "${reads_dir}/round${round}.s.fastq.gz" ]; then
                        s_size=$(stat -c%s "${reads_dir}/round${round}.s.fastq.gz" 2>/dev/null || echo 0)
                fi
                
                total_read_size=$((r1_size + r2_size + s_size))
                
                # Check for read file growth
                local growth="YES"
                if [ ${round} -gt 1 ]; then
                        if [ ${total_read_size} -eq ${prev_read_size} ]; then
                                # Only increment counter when there's truly no change
                                no_growth_cycles=$((no_growth_cycles + 1))
                                growth="NO"
                                echo "Round ${round}: No read file change detected (prev: ${prev_read_size}, curr: ${total_read_size}). No-change cycles: ${no_growth_cycles}/${max_no_growth_cycles}" >&2
                        elif [ ${total_read_size} -lt ${prev_read_size} ]; then
                                # Read set got smaller - this is good (more targeted reads)
                                no_growth_cycles=0  # Reset counter
                                growth="SMALLER"
                                echo "Round ${round}: Read file size decreased (prev: ${prev_read_size}, curr: ${total_read_size}). This is good - more targeted reads. No-change cycles reset." >&2
                        else
                                # Read set grew
                                no_growth_cycles=0  # Reset counter on growth
                                growth="YES"
                                echo "Round ${round}: Read file growth detected (prev: ${prev_read_size}, curr: ${total_read_size}). Growth cycles reset." >&2
                        fi
                else
                        echo "Round ${round}: Initial read file size: ${total_read_size} bytes" >&2
                fi
                
                # Log the size tracking information
                printf "%d\t%d\t%d\t%d\t%d\t%s\t%d\n" \
                        "${round}" "${r1_size}" "${r2_size}" "${s_size}" "${total_read_size}" "${growth}" "${no_growth_cycles}" \
                        >> "${size_tracking_log}"
                
                # NEW: Check if we should terminate due to no growth (only if we haven't found a good BLAST hit yet)
                if [ ${no_growth_cycles} -ge ${max_no_growth_cycles} ]; then
                        echo "Round ${round}: No read file growth for ${no_growth_cycles} consecutive cycles." >&2
                        echo "Read size tracking log: ${size_tracking_log}" >&2
                        # Set a flag to check after assembly/BLAST
                        local should_check_termination=true
                else
                        local should_check_termination=false
                fi
                
                prev_read_size=${total_read_size}
                
                # Attempt plasmid assembly
                echo "Round ${round}: Running SPAdes --plasmid..."
                "${SPADES_BIN}" \
                        -k ${KMER_SET} \
                        --plasmid \
                        --threads $(nproc) \
                        -1 "${reads_dir}/round${round}.1.fastq.gz" \
                        -2 "${reads_dir}/round${round}.2.fastq.gz" \
                        -s "${reads_dir}/round${round}.s.fastq.gz" \
                        -o "${outdir}" > "${log_dir}/round${round}.spades_plasmid.log" 2>&1
                
                local contigs="${outdir}/contigs.fasta"
                if [ -s "${contigs}" ] && grep -q '^>' "${contigs}"; then
                        echo "Round ${round}:    SUCCESS."
                else
                        
                        # Fallback: assemble without --plasmid and use its contigs as new reference
                        local outdir_np="${outdir}_nonplasmid"
                        mkdir -p "${outdir_np}"
                        echo "Round ${round}:    NO SUCCESS."
                        echo "Round ${round}: Running SPAdes non-plasmid..."
                        "${SPADES_BIN}" \
                                -k ${KMER_SET} \
                                --threads $(nproc) \
                                -1 "${reads_dir}/round${round}.1.fastq.gz" \
                                -2 "${reads_dir}/round${round}.2.fastq.gz" \
                                -s "${reads_dir}/round${round}.s.fastq.gz" \
                                -o "${outdir_np}" > "${log_dir}/round${round}.spades_np.log" 2>&1
                        
                        local contigs_np="${outdir_np}/contigs.fasta"
                        if ! [ -s "${contigs_np}" ] || ! grep -q '^>' "${contigs_np}"; then
                            echo "Round ${round}:    NO SUCCESS." >&2
                            
                            # NEW: Try tiered subsampling if we have a good plasmid from previous rounds
                            if [ "$found_good_np_blast_hit" = true ] && [ ${round} -gt 3 ]; then
                                echo "Round ${round}:    Assembly failed but we have a good plasmid from previous rounds. Trying tiered subsampling..." >&2
                                
                                # Try different coverage levels: 50x, 100x, 150x
                                local coverage_levels=(50 100 150)
                                local assembly_succeeded=false
                                
                                for coverage in "${coverage_levels[@]}"; do
                                    echo "Round ${round}:    Trying ${coverage}x coverage subsampling for non-plasmid assembly..." >&2
                                    
                                    # Subsample reads
                                    local subsampled_reads_dir=$(subsample_reads_for_assembly "${reads_dir}" "${round}" "${coverage}" "${sample}")
                                    
                                    # Try assembly again with subsampled reads
                                    "${SPADES_BIN}" \
                                            -k ${KMER_SET} \
                                            --threads $(nproc) \
                                            -1 "${subsampled_reads_dir}/round${round}.1.fastq.gz" \
                                            -2 "${subsampled_reads_dir}/round${round}.2.fastq.gz" \
                                            -s "${subsampled_reads_dir}/round${round}.s.fastq.gz" \
                                            -o "${outdir_np}_subsampled_${coverage}x" > "${log_dir}/round${round}.spades_np_subsampled_${coverage}x.log" 2>&1
                                    
                                    local contigs_np_sub="${outdir_np}_subsampled_${coverage}x/contigs.fasta"
                                    if [ -s "${contigs_np_sub}" ] && grep -q '^>' "${contigs_np_sub}"; then
                                        echo "Round ${round}:    ${coverage}x subsampled assembly succeeded! Using subsampled contigs." >&2
                                        contigs_np="${contigs_np_sub}"
                                        outdir_np="${outdir_np}_subsampled_${coverage}x"
                                        assembly_succeeded=true
                                        break
                                    else
                                        echo "Round ${round}:    ${coverage}x subsampled assembly failed." >&2
                                    fi
                                done
                                
                                if [ "$assembly_succeeded" != true ]; then
                                    echo "Round ${round}:    All subsampling levels failed. Sample failed." >&2
                                    mark_sample_processed "$sample" "failed"
                                    return 0
                                fi
                            else
                                echo "Round ${round}:    Assembly failed and no good plasmid from previous rounds. Sample failed." >&2
                                mark_sample_processed "$sample" "failed"
                                return 0  # Return 0 to continue with next sample instead of terminating
                            fi
                        fi
                        
                        # Perform BLAST sanity check on non-plasmid contigs
                        echo "Round ${round}:    Performing BLAST sanity check on non-plasmid contigs..."
                        local len_table_np="${outdir_np}/contig_lengths.tsv"
                        local tmp_contig_np="${outdir_np}/tmp_contig_for_blast.fasta"
                        local accepted_contig_np="${outdir_np}/accepted_contig.fasta"
                        local max_len_np=""
                        local max_name_np=""
                        local found_good_np=false
                        
                        seqkit fx2tab -nl "${contigs_np}" | sort -k2,2nr > "${len_table_np}"
                        
                        if ! [ -s "${len_table_np}" ]; then
                            echo "Round ${round}:    non-plasmid contigs length table is empty. Sample failed." >&2
                            mark_sample_processed "$sample" "failed"
                            return 0  # Return 0 to continue with next sample instead of terminating
                        fi
                        
                        while IFS=$'\t' read -r c_name c_len; do
                               # Extract this contig sequence
                               seqkit grep -n -p "${c_name}" "${contigs_np}" > "${tmp_contig_np}" 2>/dev/null || true
                              
                               if ! [ -s "${tmp_contig_np}" ]; then
                                   continue
                               fi
                              
                               # BLAST this contig against the initial reference
                               local hit_line_np
                               hit_line_np=$(blastn -query "${tmp_contig_np}" -db "${BLAST_DB_BASE}" -max_target_seqs 1 -evalue 1e-10 -outfmt '6 pident qcovs bitscore length' 2>/dev/null | head -n1 || true)
                              if [ -n "${hit_line_np}" ]; then
                                  local pid_np qcov_np bits_np len_np
                                  pid_np=$(echo "${hit_line_np}" | awk '{print $1}')
                                  qcov_np=$(echo "${hit_line_np}" | awk '{print $2}')
                                  bits_np=$(echo "${hit_line_np}" | awk '{print $3}')
                                  len_np=$(echo "${hit_line_np}" | awk '{print $4}')
                                  echo "Round ${round}:    Contig ${c_name} (len = ${c_len}) BLAST: pident = ${pid_np}%; qcovs = ${qcov_np}%; bits = ${bits_np}; alen = ${len_np}" >&2
                                  if (( $(echo "${pid_np} >= ${BLAST_PID_THRESHOLD}" | bc -l) )) && (( $(echo "${qcov_np} >= ${BLAST_QCOV_THRESHOLD}" | bc -l) )); then
                                       cp "${tmp_contig_np}" "${accepted_contig_np}"
                                       max_len_np="${c_len}"
                                       max_name_np="${c_name}"
                                       found_good_np=true
                                       found_good_blast_hit=true  # Mark that we found a good BLAST hit
                                       found_good_np_blast_hit=true  # Specifically from non-plasmid SPAdes
                                       break
                                    else
                                        echo "Round ${round}:    Contig ${c_name} failed BLAST thresholds (need >= ${BLAST_PID_THRESHOLD}% id and >= ${BLAST_QCOV_THRESHOLD}% cov)." >&2
                                    fi
                              fi
                              if [ -z "${hit_line_np}" ]; then
                                  echo "Round ${round}:    Contig ${c_name} produced no BLAST hit to reference." >&2
                              fi
                         done < "${len_table_np}"
                         
                         if [ "${found_good_np}" != true ]; then
                             # Only trigger subsampling here if a previous non-plasmid SPAdes produced a BLAST-passing contig
                             if [ "$found_good_np_blast_hit" = true ] && [ ${round} -gt 3 ]; then
                                 echo "Round ${round}:    BLAST sanity failed on non-plasmid contigs. Trying tiered subsampling..." >&2
                                 local coverage_levels=(50 100 150)
                                 local assembly_succeeded=false
                                 for coverage in "${coverage_levels[@]}"; do
                                     echo "Round ${round}:    Trying ${coverage}x coverage subsampling for non-plasmid assembly (BLAST-fail path)..." >&2
                                     local subsampled_reads_dir=$(subsample_reads_for_assembly "${reads_dir}" "${round}" "${coverage}" "${sample}")
                                     "${SPADES_BIN}" \
                                             -k ${KMER_SET} \
                                             --threads $(nproc) \
                                             -1 "${subsampled_reads_dir}/round${round}.1.fastq.gz" \
                                             -2 "${subsampled_reads_dir}/round${round}.2.fastq.gz" \
                                             -s "${subsampled_reads_dir}/round${round}.s.fastq.gz" \
                                             -o "${outdir_np}_subsampled_${coverage}x" > "${log_dir}/round${round}.spades_np_subsampled_${coverage}x.log" 2>&1
                                     local contigs_np_sub="${outdir_np}_subsampled_${coverage}x/contigs.fasta"
                                     if [ -s "${contigs_np_sub}" ] && grep -q '^>' "${contigs_np_sub}"; then
                                         echo "Round ${round}:    ${coverage}x subsampled assembly succeeded (BLAST-fail path). Using subsampled contigs." >&2
                                         contigs_np="${contigs_np_sub}"
                                         outdir_np="${outdir_np}_subsampled_${coverage}x"
                                         assembly_succeeded=true
                                         break
                                     else
                                         echo "Round ${round}:    ${coverage}x subsampled assembly failed (BLAST-fail path)." >&2
                                     fi
                                 done
                                 if [ "${assembly_succeeded}" != true ]; then
                                     echo "Round ${round}:    All subsampling levels failed after BLAST sanity failure. Sample failed." >&2
                                     mark_sample_processed "$sample" "failed"
                                     return 0
                                 fi
                                 # If subsampled assembly succeeded, redo BLAST sanity on new contigs
                                 echo "Round ${round}:    Re-running BLAST sanity on subsampled non-plasmid contigs..." >&2
                                 len_table_np="${outdir_np}/contig_lengths.tsv"
                                 tmp_contig_np="${outdir_np}/tmp_contig_for_blast.fasta"
                                 accepted_contig_np="${outdir_np}/accepted_contig.fasta"
                                 found_good_np=false
                                 seqkit fx2tab -nl "${contigs_np}" | sort -k2,2nr > "${len_table_np}"
                                 echo "Round ${round}:    BLAST thresholds: pident >= ${BLAST_PID_THRESHOLD} and qcovs >= ${BLAST_QCOV_THRESHOLD} (subsampled)" >&2
                                 while IFS=$'\t' read -r c_name c_len; do
                                     seqkit grep -n -p "${c_name}" "${contigs_np}" > "${tmp_contig_np}" 2>/dev/null || true
                                     if ! [ -s "${tmp_contig_np}" ]; then
                                         continue
                                     fi
                                     local hit_line_np
                                     hit_line_np=$(blastn -query "${tmp_contig_np}" -db "${BLAST_DB_BASE}" -max_target_seqs 1 -evalue 1e-10 -outfmt '6 pident qcovs bitscore length' 2>/dev/null | head -n1 || true)
                                     if [ -n "${hit_line_np}" ]; then
                                         local pid_np qcov_np bits_np len_np
                                         pid_np=$(echo "${hit_line_np}" | awk '{print $1}')
                                         qcov_np=$(echo "${hit_line_np}" | awk '{print $2}')
                                         bits_np=$(echo "${hit_line_np}" | awk '{print $3}')
                                         len_np=$(echo "${hit_line_np}" | awk '{print $4}')
                                         echo "Round ${round}:    NP-sub contig ${c_name} (len = ${c_len}) BLAST: pident = ${pid_np}%; qcovs = ${qcov_np}%; bits = ${bits_np}; alen = ${len_np}" >&2
                                         if (( $(echo "${pid_np} >= ${BLAST_PID_THRESHOLD}" | bc -l) )) && (( $(echo "${qcov_np} >= ${BLAST_QCOV_THRESHOLD}" | bc -l) )); then
                                             cp "${tmp_contig_np}" "${accepted_contig_np}"
                                             max_len_np="${c_len}"
                                             max_name_np="${c_name}"
                                             found_good_np=true
                                             found_good_blast_hit=true
                                             found_good_np_blast_hit=true  # keep consistency: success from non-plasmid SPAdes
                                             break
                                         else
                                             echo "Round ${round}:    NP-sub contig ${c_name} failed BLAST thresholds (need >= ${BLAST_PID_THRESHOLD}% id and >= ${BLAST_QCOV_THRESHOLD}% cov)." >&2
                                         fi
                                     fi
                                     if [ -z "${hit_line_np}" ]; then
                                         echo "Round ${round}:    NP-sub contig ${c_name} produced no BLAST hit to reference." >&2
                                     fi
                                 done < "${len_table_np}"
                                 if [ "${found_good_np}" != true ]; then
                                     echo "Round ${round}:    Subsampled non-plasmid contigs still failed BLAST sanity. Sample failed." >&2
                                     mark_sample_processed "$sample" "failed"
                                     return 0
                                 fi
                             else
                                 echo "Round ${round}:    No non-plasmid contigs passed BLAST sanity check (>=${BLAST_PID_THRESHOLD}% identity and >=${BLAST_QCOV_THRESHOLD}% coverage). Sample failed." >&2
                                 mark_sample_processed "$sample" "failed"
                                 return 0
                             fi
                         fi
                         
                         # Use the accepted non-plasmid contig as the new reference
                          current_ref="${accepted_contig_np}"
                          
                          # Add to plasmid size history
                          plasmid_size_history+=("${max_len_np}")
                          
                          # Check if we should terminate due to no growth
                          if [ "${should_check_termination}" = true ]; then
                                  if [ "${found_good_blast_hit}" = true ]; then
                                          echo "Round ${round}: No read file growth detected, but good BLAST hit found. Using current contig as final plasmid." >&2
                                          largest_fa="${accepted_contig_np}"
                                          prev_max_len=${max_len_np}
                                          break
                                  else
                                          echo "Round ${round}: No read file growth for ${no_growth_cycles} consecutive cycles AND no good BLAST hits found. Terminating sample as failed." >&2
                                          mark_sample_processed "$sample" "failed"
                                          return 0
                                  fi
                          fi
                          
                          # NEW: Check for plateau in plasmid size
                          if detect_plateau "${plasmid_size_history[@]}"; then
                                  echo "Round ${round}: Plateau detected in plasmid size growth. Terminating optimization." >&2
                                  largest_fa="${accepted_contig_np}"
                                  prev_max_len=${max_len_np}
                                  break
                          fi
                          
                          round=$((round + 1))
                          continue
                fi
                
                # Compute contig lengths and perform BLAST sanity check over descending lengths
		local len_table="${outdir}/contig_lengths.tsv"
		local tmp_contig_fa="${outdir}/tmp_contig_for_blast.fasta"
		local accepted_contig_fa="${outdir}/accepted_contig.fasta"
		local max_len=""
		local max_name=""
		local found_good=false
		
		seqkit fx2tab -nl "${contigs}" | sort -k2,2nr > "${len_table}"
		
		if ! [ -s "${len_table}" ]; then
			echo "Round ${round}: contigs length table is empty. Sample failed." >&2
			mark_sample_processed "$sample" "failed"
			return 0  # Return 0 to continue with next sample instead of terminating
		fi
		
		echo "Round ${round}: BLAST thresholds: pident >= ${BLAST_PID_THRESHOLD} and qcovs >= ${BLAST_QCOV_THRESHOLD} (plasmid path)" >&2
		while IFS=$'\t' read -r c_name c_len; do
			# Extract this contig sequence
			
			seqkit grep -n -p "${c_name}" "${contigs}" > "${tmp_contig_fa}" 2>/dev/null || true
			
			if ! [ -s "${tmp_contig_fa}" ]; then
				continue
			fi
			
						# BLAST this contig against the initial reference
			local hit_line
			hit_line=$(blastn -query "${tmp_contig_fa}" -db "${BLAST_DB_BASE}" -max_target_seqs 1 -evalue 1e-10 -outfmt '6 pident qcovs bitscore length' 2>/dev/null | head -n1 || true)
			if [ -n "${hit_line}" ]; then
				local pid qcov bits len
				pid=$(echo "${hit_line}" | awk '{print $1}')
				qcov=$(echo "${hit_line}" | awk '{print $2}')
				bits=$(echo "${hit_line}" | awk '{print $3}')
				len=$(echo "${hit_line}" | awk '{print $4}')
								if (( $(echo "${pid} >= ${BLAST_PID_THRESHOLD}" | bc -l) )) && (( $(echo "${qcov} >= ${BLAST_QCOV_THRESHOLD}" | bc -l) )); then
 					cp "${tmp_contig_fa}" "${accepted_contig_fa}"
 					max_len="${c_len}"
 					max_name="${c_name}"
 					found_good=true
 					found_good_blast_hit=true  # Mark that we found a good BLAST hit
 					break
 				fi
 			fi
		done < "${len_table}"
		
		if [ "${found_good}" != true ]; then
			echo "Round ${round}: No contigs passed BLAST sanity check (>=${BLAST_PID_THRESHOLD}% identity and >=${BLAST_QCOV_THRESHOLD}% coverage). Sample failed." >&2
			mark_sample_processed "$sample" "failed"
			return 0  # Return 0 to continue with next sample instead of terminating
		fi
		
		# Use the accepted contig as the largest/reference for next round
		largest_fa="${accepted_contig_fa}"
		printf "%s\n" "${max_len}" > "${outdir}/max_contig_len.txt"
		echo "Round ${round}: selected contig length = ${max_len} (name: ${max_name}) after BLAST sanity check."
		
		if ! [ -s "${largest_fa}" ]; then
			echo "Round ${round}: failed to extract accepted contig sequence. Sample failed." >&2
			mark_sample_processed "$sample" "failed"
			return 0  # Return 0 to continue with next sample instead of terminating
		fi
                
                # Add to plasmid size history
                plasmid_size_history+=("${max_len}")
                
                # Update prev_max_len with current max_len before checking plateau
                prev_max_len=${max_len}
                
                # Check if we should terminate due to no growth
                if [ "${should_check_termination}" = true ]; then
                        if [ "${found_good_blast_hit}" = true ]; then
                                echo "Round ${round}: No read file growth detected, but good BLAST hit found. Using current contig as final plasmid." >&2
                                break
                        else
                                echo "Round ${round}: No read file growth for ${no_growth_cycles} consecutive cycles AND no good BLAST hits found. Terminating sample as failed." >&2
                                mark_sample_processed "$sample" "failed"
                                return 0
                        fi
                fi
                
                # NEW: Check for plateau in plasmid size (replaces max rounds logic)
                if detect_plateau "${plasmid_size_history[@]}"; then
                        echo "Round ${round}: Plateau detected in plasmid size growth. Terminating optimization." >&2
                        plateau_detected=true
                        break
                fi
                
                # REMOVED: Max rounds check - now we rely on plateau detection
                # The script will continue until plateau is detected or max rounds is reached as a safety net
                
                current_ref="${largest_fa}"
                round=$((round + 1))
        done
        
        echo "============================================================================================================" >&2
        if [ ${round} -gt ${MAX_ROUNDS} ]; then
                echo "========= Sample ${sample} complete (MAX ROUNDS REACHED). Final largest contig length: ${prev_max_len} =========" >&2
        elif [ "${plateau_detected}" = true ]; then
                echo "========= Sample ${sample} complete (PLATEAU DETECTED). Final largest contig length: ${prev_max_len} =========" >&2
        elif [ ${no_growth_cycles} -ge ${max_no_growth_cycles} ] && [ "${found_good_blast_hit}" = true ]; then
                echo "========= Sample ${sample} complete (NO GROWTH BUT GOOD BLAST HIT). Final largest contig length: ${prev_max_len} =========" >&2
        elif [ ${no_growth_cycles} -ge ${max_no_growth_cycles} ]; then
                echo "========= Sample ${sample} complete (NO GROWTH DETECTED). Terminated after ${no_growth_cycles} no-growth cycles =========" >&2
        else
                echo "========= Sample ${sample} complete (UNKNOWN REASON). Final largest contig length: ${prev_max_len} =========" >&2
        fi
        echo "============================================================================================================" >&2
        
        # Copy the best plasmid contig to the sample directory for easy access
        if [ -n "${largest_fa}" ] && [ -s "${largest_fa}" ]; then
                # Copy with renamed header to plasmids directory
                sed "s/^>.*/>${sample}/" "${largest_fa}" > "${PLASMIDS_DIR}/${sample}_plasmid.fasta"
                echo "Renamed plasmid copied to: ${PLASMIDS_DIR}/${sample}_plasmid.fasta" >&2
                
                # Copy the last used reads to the reads directory
                # Prefer current round if files exist (loop exited via break), else fall back to previous round
                local final_round=${round}
                if ! [ -f "${reads_dir}/round${final_round}.1.fastq.gz" ] || ! [ -f "${reads_dir}/round${final_round}.2.fastq.gz" ] || ! [ -f "${reads_dir}/round${final_round}.s.fastq.gz" ]; then
                        final_round=$((round - 1))
                fi
                if [ ${final_round} -gt 0 ] && \
                   [ -f "${reads_dir}/round${final_round}.1.fastq.gz" ] && \
                   [ -f "${reads_dir}/round${final_round}.2.fastq.gz" ] && \
                   [ -f "${reads_dir}/round${final_round}.s.fastq.gz" ]; then
                        cp "${reads_dir}/round${final_round}.1.fastq.gz" "${READS_DIR}/${sample}_final_round_1.fastq.gz"
                        cp "${reads_dir}/round${final_round}.2.fastq.gz" "${READS_DIR}/${sample}_final_round_2.fastq.gz"
                        cp "${reads_dir}/round${final_round}.s.fastq.gz" "${READS_DIR}/${sample}_final_round_s.fastq.gz"
                        echo "Final round reads copied to reads directory from round ${final_round}:" >&2
                        echo "  - ${READS_DIR}/${sample}_final_round_1.fastq.gz" >&2
                        echo "  - ${READS_DIR}/${sample}_final_round_2.fastq.gz" >&2
                        echo "  - ${READS_DIR}/${sample}_final_round_s.fastq.gz" >&2
                else
                        echo "WARNING: Unable to locate final round read files for copying (checked rounds ${round} and $((round-1)))." >&2
                fi
                
                # Mark sample as successfully processed
                mark_sample_processed "$sample" "success"
        else
                # Mark sample as failed
                mark_sample_processed "$sample" "failed"
        fi
        
        return 0  # Always return success to continue with next sample
}

# Process all samples
total_samples=${#BATCH_SAMPLES[@]}
current_sample_num=1
failed_samples=0
skipped_samples=0

# Check required tools before starting processing
echo "Checking required tools..." >&2
if ! check_required_tools; then
        exit 1
fi

echo "======================================================" >&2
echo "PROCESSING STATUS:" >&2
echo "======================================================" >&2

for sample in "${BATCH_SAMPLES[@]}"; do
        echo "======================================================" >&2
        echo "=== Processing sample ${current_sample_num} of ${total_samples}" >&2
        
        # Check if sample was already processed
        if is_sample_processed "$sample"; then
                echo "Sample ${sample} was already processed. Skipping." >&2
                skipped_samples=$((skipped_samples + 1))
        else
                # Process the sample and capture the exit code
                echo "Starting to process sample: ${sample}" >&2
                # Use set +e to temporarily disable exit on error for this function call
                set +e
                if process_sample "$sample"; then
                        # Sample processed successfully
                        echo "Sample ${sample} processed successfully." >&2
                else
                        echo "Failed to process sample: ${sample}" >&2
                        failed_samples=$((failed_samples + 1))
                fi
                set -e  # Re-enable exit on error
        fi
        
        current_sample_num=$((current_sample_num + 1))
        echo "Completed sample ${current_sample_num} of ${total_samples}" >&2
done

echo "======================================================" >&2
echo "BATCH PROCESSING COMPLETE!" >&2
echo "======================================================" >&2
echo "Total samples: ${total_samples}" >&2
echo "Skipped samples (already processed): ${skipped_samples}" >&2
echo "Failed samples: ${failed_samples}" >&2
echo "Successfully processed: $((total_samples - skipped_samples - failed_samples))" >&2

echo "" >&2
echo "============================================================================================================" >&2
echo "RESULTS SUMMARY:" >&2
echo "============================================================================================================" >&2
echo "Top-level results directory: ${RESULTS_DIR}/" >&2
echo "Working directory (individual samples): ${WORKING_DIR}/" >&2
echo "Plasmids directory: ${PLASMIDS_DIR}/" >&2
echo "Reads directory: ${READS_DIR}/" >&2
echo "" >&2
echo "Final plasmids with renamed headers: ${PLASMIDS_DIR}/*_plasmid.fasta" >&2
echo "Final round reads: ${READS_DIR}/*_final_round_*.fastq.gz" >&2
echo "Sample-specific working directories: ${WORKING_DIR}/{SAMPLE}/" >&2
echo "" >&2
echo "Processing status file: ${PROCESSING_STATUS_FILE}" >&2
if [ -s "${FAILED_SAMPLES_FILE}" ]; then
        echo "Failed samples: ${FAILED_SAMPLES_FILE}" >&2
fi
if [ -s "${SUCCESS_SAMPLES_FILE}" ]; then
        echo "Successfully processed samples: ${SUCCESS_SAMPLES_FILE}" >&2
fi
echo "============================================================================================================" >&2

# Exit with error code if any samples failed
if [ ${failed_samples} -gt 0 ]; then
        exit 1
fi