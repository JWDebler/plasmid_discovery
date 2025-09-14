#!/bin/bash

# VLE Mapper Script
# This script downloads SRA entries, maps reads to a reference genome,
# and identifies entries with sufficient coverage (>1x) as potential VLEs.
# 
# New functionality: During read mapping, only reads mapping to the template
# are retained, and alignment files (BAM) are saved to reference-named
# folders for successful mappings with >1x coverage.
# 
# CSV processing enhancements: Supports SampleName, LibraryName,
# and CenterName fields. Sample names and library names are preserved
# exactly as they appear in the CSV to maintain data integrity
# (e.g., "AR0037_AR0037" stays intact).

# Debug mode flag
DEBUG_MODE=false

# Color codes
RED='\033[0;31m'
NC='\033[0m' # No Color

# Process tracking arrays
declare -A SUBPROCESS_PIDS
declare -A SUBPROCESS_NAMES

# Function to add subprocess to tracking
track_subprocess() {
    local pid="$1"
    local name="$2"
    SUBPROCESS_PIDS[$pid]=$pid
    SUBPROCESS_NAMES[$pid]=$name
    debug_log "Started subprocess: $name (PID: $pid)"
}

# Function to remove subprocess from tracking
untrack_subprocess() {
    local pid="$1"
    if [[ -n "${SUBPROCESS_PIDS[$pid]}" ]]; then
        debug_log "Completed subprocess: ${SUBPROCESS_NAMES[$pid]} (PID: $pid)"
        unset SUBPROCESS_PIDS[$pid]
        unset SUBPROCESS_NAMES[$pid]
    fi
}

# Function to cleanup all subprocesses
cleanup_subprocesses() {
    info_log "Cleaning up subprocesses..."
    
    # First, try graceful termination
    for pid in "${!SUBPROCESS_PIDS[@]}"; do
        if kill -0 "$pid" 2>/dev/null; then
            debug_log "Terminating ${SUBPROCESS_NAMES[$pid]} (PID: $pid)"
            kill -TERM "$pid" 2>/dev/null
        fi
    done
    
    # Wait a moment for processes to terminate
    sleep 1
    
    # Force kill any remaining processes
    for pid in "${!SUBPROCESS_PIDS[@]}"; do
        if kill -0 "$pid" 2>/dev/null; then
            debug_log "Force killing ${SUBPROCESS_NAMES[$pid]} (PID: $pid)"
            kill -9 "$pid" 2>/dev/null
        fi
    done
    
    # Clear the tracking arrays
    SUBPROCESS_PIDS=()
    SUBPROCESS_NAMES=()
    
    # Kill any remaining fasterq-dump or mapping processes
    pkill -9 -f "fasterq-dump" 2>/dev/null
    pkill -9 -f "bwa-mem2" 2>/dev/null
    pkill -9 -f "minimap2" 2>/dev/null
    pkill -9 -f "samtools" 2>/dev/null
}

# Function to print debug messages
debug_log() {
    if [[ "$DEBUG_MODE" == "true" ]]; then
        echo -n "[DEBUG] $*" >&2
        echo >&2
    fi
}

# Function to print info messages
info_log() {
    echo "[INFO]  $*"
}

# Function to print error messages
error_log() {
    echo -e "${RED}[ERROR] $*${NC}" >&2
}

# Function to print usage information
print_usage() {
    cat << EOF
Usage: $(basename "$0") [options]

Options:
    -h, --help          Show this help message and exit
    -c, --csv FILE      Path to input CSV file (optional)
    -d, --db FILE       Path to SQLite database file (optional)
                        If not provided, creates new database as {organism}sra_wgs.db
                        If provided with CSV, updates non-pending entries
    -o, --organism STR  Organism to search for in SRA (default: fungi)
                        Used only when fetching new data from SRA
    -x, --connections INT  Number of concurrent connections to SRA (default: 5)
    -r, --retry         Retry failed downloads only (selective retry)
                        Only retries entries that failed during download/extraction phase
    -R, --retry-all     Retry ALL failed entries (comprehensive retry)
                        Resets all failed entries regardless of failure reason
    -f, --reference FILE Path to reference genome FASTA file (required)
                        Downloaded SRA entries will be extracted, mapped to the reference,
                        filtered for mapped reads only, and analyzed for coverage.
                        Alignment files (BAM) are saved to 'alignments_mapping_to_reference'
                        folder for entries with >1x coverage.
                        CSV processing supports SampleName, LibraryName, and CenterName fields.
                        Sample names and library names are preserved exactly as they appear
                        in the CSV (including underscores).
    -D, --debug        Enable debug mode for verbose output

Examples:
    # Create new database from SRA query (default: fungi)
    $(basename "$0") -f reference.fasta
    
    # Create new database for specific organism
    $(basename "$0") -o "Aspergillus" -f reference.fasta
    
    # Create new database from CSV file
    $(basename "$0") -c existing.csv -f reference.fasta
    
    # Resume processing existing database
    $(basename "$0") -d existing.db -f reference.fasta
    
    # Resume and retry only download failures
    $(basename "$0") -d existing.db -f reference.fasta -r
    
    # Resume and retry ALL failed entries (recommended for systematic fixes)
    $(basename "$0") -d existing.db -f reference.fasta -R
    
    # Update existing database with new entries from CSV file
    $(basename "$0") -c new_entries.csv -d existing.db -f reference.fasta
    
    # Update existing database with new entries for a specific organism
    $(basename "$0") -d existing.db -o "Aspergillus" -f reference.fasta
    
    # Run with debug output enabled
    $(basename "$0") -d existing.db -f reference.fasta -D

Retry Options Explained:
    -r (--retry):        Selective retry - only retries entries that failed during
                        download/prefetch/extraction. Use when you had network issues
                        or SRA server problems.
    
    -R (--retry-all):    Comprehensive retry - retries ALL failed entries regardless
                        of why they failed. Use when you've fixed a systematic issue
                        (missing tools, wrong reference file, etc.) or want to give
                        everything a fresh start.

Output Directories:
    - {organism}/                     # Top-level folder named after the organism from -o option
        - {species}/                  # Subfolder named after each species with matches
            - reads/                  # Compressed FASTQ files for entries with >1x coverage
            - alignments/             # BAM alignment files for mapped reads with >1x coverage
    - {organism}sra_wgs.db            # SQLite database containing all SRA entries and status
EOF
}

# Function to classify sequencing technology
classify_technology() {
    local platform="$1"
    local model="$2"
    
    # Convert to uppercase for consistent matching
    platform=$(echo "$platform" | tr '[:lower:]' '[:upper:]')
    model=$(echo "$model" | tr '[:lower:]' '[:upper:]')
    
    # Illumina and similar short-read platforms
    if [[ "$platform" == "ILLUMINA" ]] || 
       [[ "$platform" == "DNBSEQ" ]] || 
       [[ "$platform" =~ ^BGI ]] || 
       [[ "$platform" == "MGISEQ" ]] || 
       [[ "$platform" == "BGISEQ" ]] || 
       [[ "$platform" == "COMPLETE_GENOMICS" ]] ||  
       [[ "$model" == *"ILLUMINA"* ]] ||
       [[ "$model" == *"HISEQ"* ]] ||
       [[ "$model" == *"NEXTSEQ"* ]] ||
       [[ "$model" == *"NOVASEQ"* ]] ||
       [[ "$model" == *"MISEQ"* ]] ||
       [[ "$model" == *"ISEQ"* ]]; then
        echo "illumina"
        return 0
    fi
    
    # Oxford Nanopore platforms
    if [[ "$platform" == "OXFORD_NANOPORE" ]] || 
       [[ "$platform" == "ONT" ]] || 
       [[ "$model" =~ (MINION|GRIDION|PROMETHION) ]]; then
        echo "nanopore"
        return 0
    fi
    
    # PacBio platforms
    if [[ "$platform" == "PACBIO" ]] || 
       [[ "$platform" =~ ^PACBIO ]] || 
       [[ "$model" =~ (SEQUEL|REVIO|RS II) ]]; then
        echo "pacbio"
        return 0
    fi
    
    # Ion Torrent platforms - moved up for clarity and priority
    if [[ "$platform" == "ION_TORRENT" ]] || 
       [[ "$model" =~ (ION|PROTON|PGM|S5|GENEXUS) ]]; then
        echo "iontorrent"
        return 0
    fi
    
    # 454 pyrosequencing platforms - moved up for clarity and priority
    if [[ "$platform" == "LS454" ]] || 
       [[ "$model" =~ (454|GS FLX|GS JUNIOR) ]]; then
        echo "454"
        return 0
    fi
    
    # Capillary sequencing platforms
    if [[ "$platform" == "CAPILLARY" ]] || 
       [[ "$model" =~ (3730XL|3730|3130XL) ]]; then
        echo "capillary"
        return 0
    fi
    
    # ABI SOLiD platforms
    if [[ "$platform" == "ABI_SOLID" ]] || 
       [[ "$model" =~ (SOLID|5500XL) ]]; then
        echo "solid"
        return 0
    fi
    
    # Unknown/other technology
    echo "other"
    return 1
}

# Function to sanitize string for SQL
sanitize_sql_string() {
    local str="$1"
    # Replace single quotes with two single quotes (SQL escape)
    echo "${str//\'/\'\'}"
}

# Function to sanitize filename
sanitize_filename() {
    local name="$1"
    # Convert to lowercase, replace spaces and special chars with underscore, then compress multiple underscores
    echo "$name" | tr '[:upper:]' '[:lower:]' | tr -c '[:alnum:]' '_' | tr -s '_' | sed 's/^_\|_$//g'
}

# Function to execute SQLite command with retries
execute_sqlite_cmd() {
    local db_file="$1"
    local cmd="$2"
    local max_retries=5
    local retry_delay=1
    local attempt=1
    local lock_file="${db_file}.lock"
    local result=""
    local suppress_debug="${3:-false}"  # New parameter to suppress debug output
    
    # Ensure the database directory exists
    local db_dir=$(dirname "$db_file")
    if [[ ! -d "$db_dir" ]]; then
        mkdir -p "$db_dir" || {
            error_log "Failed to create database directory: $db_dir"
            return 1
        }
    fi
    
    # Check if we can write to the database directory
    if [[ ! -w "$db_dir" ]]; then
        error_log "Database directory is not writable: $db_dir"
        return 1
    fi
    
    # If database doesn't exist, try to create it
    if [[ ! -f "$db_file" ]]; then
        debug_log "Creating new database file: $db_file"
        sqlite3 "$db_file" "SELECT 1;" >/dev/null 2>/dev/null || {
            error_log "Failed to create database file: $db_file"
            return 1
        }
    fi
    
    # Check if database is writable
    if [[ -f "$db_file" ]] && [[ ! -w "$db_file" ]]; then
        error_log "Database file is not writable: $db_file"
        return 1
    fi
    
    while ((attempt <= max_retries)); do
        # Use a temporary file for output
        local tmp_output=$(mktemp)
        local tmp_error=$(mktemp)
        
        # Only log if we're not suppressing debug output
        if [[ "$suppress_debug" != "true" ]]; then
            debug_log "Executing SQLite command (attempt $attempt/$max_retries)..."
        fi
        
        if (
            flock -x -w 10 200 || exit 1
            sqlite3 "$db_file" "$cmd" > "$tmp_output" 2>"$tmp_error"
        ) 200>"$lock_file"; then
            # Command succeeded
            result=$(<"$tmp_output")
            rm -f "$tmp_output" "$tmp_error" "$lock_file"
            echo "$result"
            return 0
        else
            local error_msg=$(<"$tmp_error")
            
            # Only log errors if we're not suppressing debug output
            if [[ "$suppress_debug" != "true" ]]; then
                debug_log "SQLite error (attempt $attempt): $error_msg"
            
                # Check for specific error conditions
                if [[ "$error_msg" == *"database is locked"* ]]; then
                    debug_log "Database is locked, waiting before retry..."
                elif [[ "$error_msg" == *"disk I/O error"* ]]; then
                    error_log "Disk I/O error occurred, check disk space and permissions"
                    rm -f "$tmp_output" "$tmp_error" "$lock_file"
                    return 1
                elif [[ "$error_msg" == *"permission denied"* ]]; then
                    error_log "Permission denied accessing database"
                    rm -f "$tmp_output" "$tmp_error" "$lock_file"
                    return 1
                fi
            else
                # For critical errors, always log even if suppressing debug output
                if [[ "$error_msg" == *"disk I/O error"* ]]; then
                    error_log "Disk I/O error occurred, check disk space and permissions"
                    rm -f "$tmp_output" "$tmp_error" "$lock_file"
                    return 1
                elif [[ "$error_msg" == *"permission denied"* ]]; then
                    error_log "Permission denied accessing database"
                    rm -f "$tmp_output" "$tmp_error" "$lock_file"
                    return 1
                fi
            fi
            
            rm -f "$tmp_output" "$tmp_error"
            
            if ((attempt == max_retries)); then
                error_log "Failed to execute SQLite command after $max_retries attempts"
                error_log "Last error: $error_msg"
                rm -f "$lock_file"
                return 1
            fi
            
            sleep $retry_delay
            ((retry_delay *= 2))
        fi
        ((attempt++))
    done
    
    rm -f "$lock_file"
    return 1
}

# Function to initialize database
init_database() {
    local db_file="$1"
    info_log "Initializing database $db_file..."
    
    # Create the database if it doesn't exist
    if ! execute_sqlite_cmd "$db_file" "
        CREATE TABLE IF NOT EXISTS sra_entries (
            sra_id TEXT PRIMARY KEY,
            sample_name TEXT,
            library_name TEXT,
            organism TEXT,
            platform TEXT,
            model TEXT,
            technology TEXT,
            bases INTEGER,
            center_name TEXT,
            message TEXT,
            status TEXT DEFAULT 'pending'
        );"; then
        error_log "Failed to initialize database schema"
        return 1
    fi
    
    # Verify the table was created
    if ! execute_sqlite_cmd "$db_file" "SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='sra_entries';" | grep -q "1"; then
        error_log "Failed to verify database initialization"
        return 1
    fi
    
    debug_log "Database initialized successfully"
    return 0
}

# Function to parse CSV line with proper handling of quoted fields
parse_csv_line() {
    local line="$1"
    local -n result="$2"  # Pass array by reference
    
    local field=""
    local in_quotes=false
    local i=0
    
    # Process the line character by character
    for (( pos=0; pos<${#line}; pos++ )); do
        local char="${line:$pos:1}"
        
        if [[ "$char" == '"' ]]; then
            # Toggle quote state
            in_quotes=$([ "$in_quotes" = true ] && echo false || echo true)
            field+="$char"  # Preserve quotes for now
        elif [[ "$char" == ',' && "$in_quotes" == false ]]; then
            # This is a field separator
            result[$i]="$field"
            field=""
            ((i++))
        else
            # Normal character
            field+="$char"
        fi
    done
    
    # Add the last field
    result[$i]="$field"
    
    # Clean up quoted fields by removing enclosing quotes
    for (( j=0; j<=i; j++ )); do
        if [[ "${result[$j]}" =~ ^\".*\"$ ]]; then
            result[$j]="${result[$j]:1:${#result[$j]}-2}"
        fi
    done
    
    return 0
}

# Function to process CSV file
process_csv() {
    local csv_file="$1"
    local db_file="$2"
    local processed=0
    local skipped=0
    local illumina_count=0
    local nanopore_count=0
    local pacbio_count=0
    local iontorrent_count=0
    local ls454_count=0
    local capillary_count=0
    local solid_count=0
    local other_count=0
    local batch_size=100
    local current_batch=0
    local temp_sql="temp_commands.sql"
    
    info_log "Processing CSV file: $csv_file"
    
    # Process header to get column positions
    local header
    read -r header < "$csv_file"
    declare -a headers
    parse_csv_line "$header" headers
    
    # Find column indices
    local run_idx=-1 sample_idx=-1 library_idx=-1 org_idx=-1 platform_idx=-1 model_idx=-1 bases_idx=-1 center_idx=-1
    for i in "${!headers[@]}"; do
        case "${headers[$i]}" in
            "Run") run_idx=$i ;;
            "SampleName") sample_idx=$i ;;
            "LibraryName") library_idx=$i ;;
            "ScientificName") org_idx=$i ;;
            "Platform") platform_idx=$i ;;
            "Model") model_idx=$i ;;
            "bases") bases_idx=$i ;;
            "CenterName") center_idx=$i ;;
        esac
    done
    
    # Verify all required columns were found
    local missing_columns=()
    if [[ $run_idx -eq -1 ]]; then missing_columns+=("Run"); fi
    if [[ $sample_idx -eq -1 ]]; then missing_columns+=("SampleName"); fi
    if [[ $library_idx -eq -1 ]]; then missing_columns+=("LibraryName"); fi
    if [[ $org_idx -eq -1 ]]; then missing_columns+=("ScientificName"); fi
    if [[ $platform_idx -eq -1 ]]; then missing_columns+=("Platform"); fi
    if [[ $model_idx -eq -1 ]]; then missing_columns+=("Model"); fi
    if [[ $bases_idx -eq -1 ]]; then missing_columns+=("bases"); fi
    
    if [[ ${#missing_columns[@]} -gt 0 ]]; then
        error_log "Error: Missing required columns in CSV file: ${missing_columns[*]}"
        error_log "Available columns: ${headers[*]}"
        return 1
    fi
    
    # Start the first transaction
    {
        echo "BEGIN TRANSACTION;"
    } > "$temp_sql"
    
    # Process each line in batches
    while IFS= read -r line; do
        # Skip header
        if ((processed == 0)); then
            ((processed++))
            continue
        fi
        
        # Parse CSV line with proper handling of quoted fields
        declare -a fields
        parse_csv_line "$line" fields
        
        # Extract and sanitize fields
        local run="${fields[$run_idx]}"
        local sample="${fields[$sample_idx]}"
        local library="${fields[$library_idx]}"
        local organism="${fields[$org_idx]}"
        local platform="${fields[$platform_idx]}"
        local model="${fields[$model_idx]}"
        local bases="${fields[$bases_idx]}"
        local center="${fields[$center_idx]:-}"  # Default to empty string if not found
        
        # Skip if any required field is empty or bases is 0
        if [[ -z "$run" ]] || [[ -z "$platform" ]] || [[ -z "$bases" ]] || ((bases == 0)); then
            ((skipped++))
            continue
        fi
        
        # Classify technology
        local tech
        tech=$(classify_technology "$platform" "$model")
        
        # Update counts
        case "$tech" in
            "illumina") ((illumina_count++)) ;;
            "nanopore") ((nanopore_count++)) ;;
            "pacbio") ((pacbio_count++)) ;;
            "iontorrent") ((iontorrent_count++)) ;;
            "454") ((ls454_count++)) ;;
            "capillary") ((capillary_count++)) ;;
            "solid") ((solid_count++)) ;;
            "other") ((other_count++)) ;;
        esac
        
        # Sanitize inputs (but preserve sample name as-is to keep underscores)
        run=$(sanitize_sql_string "$run")
        # Don't sanitize sample name - preserve it exactly as it appears in CSV
        sample=$(sanitize_sql_string "$sample")
        organism=$(sanitize_sql_string "$organism")
        platform=$(sanitize_sql_string "$platform")
        model=$(sanitize_sql_string "$model")
        center=$(sanitize_sql_string "$center")
        
        # Append INSERT statement to temp file
        {
            echo "INSERT OR IGNORE INTO sra_entries (sra_id, sample_name, library_name, organism, platform, model, technology, bases, center_name, status) 
                  VALUES ('$run', '$sample', '$library', '$organism', '$platform', '$model', '$tech', $bases, '$center', 'pending');"
        } >> "$temp_sql"
        
        ((current_batch++))
        
        # Process batch if we've reached batch_size
        if ((current_batch >= batch_size)); then
            # Commit the current batch
            {
                echo "COMMIT;"
            } >> "$temp_sql"
            
            # Execute current batch - suppress debug output and redirect all output to /dev/null
            if ! execute_sqlite_cmd "$db_file" "$(cat $temp_sql)" "true" > /dev/null 2>&1; then
                error_log "Failed to process batch ending at line $processed"
                rm -f "$temp_sql"
                return 1
            fi
            
            # Start new batch
            {
                echo "BEGIN TRANSACTION;"
            } > "$temp_sql"
            
            current_batch=0
        fi
        
        ((processed++))
        
        # Show progress
        if ((processed % 1000 == 0)); then
            info_log "Processed $processed entries (Skipped: $skipped)..."
        fi
    done < "$csv_file"
    
    # Process final batch if any
    if ((current_batch > 0)); then
        {
            echo "COMMIT;"
        } >> "$temp_sql"
        
        if ! execute_sqlite_cmd "$db_file" "$(cat $temp_sql)" "true" > /dev/null 2>&1; then
            error_log "Failed to process final batch"
            rm -f "$temp_sql"
            return 1
        fi
    fi
    
    rm -f "$temp_sql"
    
    # Verify the entries were actually inserted - don't suppress debug for the count query
    local actual_count
    actual_count=$(execute_sqlite_cmd "$db_file" "SELECT COUNT(*) FROM sra_entries;" 2>/dev/null)
    
    if [[ -z "$actual_count" ]] || [[ "$actual_count" -eq 0 ]]; then
        error_log "Failed to insert entries into the database"
        return 1
    fi
    
    # Print summary
    info_log "Processing complete!"
    info_log "Total entries processed: $((processed-1))"
    info_log "Entries skipped: $skipped"
    info_log "Technology distribution:"
    info_log "Illumina: $illumina_count"
    info_log "Nanopore: $nanopore_count"
    info_log "PacBio: $pacbio_count"
    info_log "IonTorrent: $iontorrent_count"
    info_log "454: $ls454_count"
    info_log "Capillary: $capillary_count"
    info_log "SOLiD: $solid_count"
    info_log "Other: $other_count"
    
    return 0
}

# Function to process SRA entry (download, extract, and map)
process_sra_entry() {
    local sra_id="$1"
    local db_file="$2"
    local slot="$3"
    local counter_file="$4"
    local total_entries="$5"
    local reference="$6"
    local csv_organism_name="$7"  # New parameter for organism name from CSV filename
    local threads="$(nproc)"
    
    # Extract reference filename without path and extension for folder naming
    local ref_basename=$(basename "$reference")
    local ref_name="${ref_basename%.*}"  # Remove extension
    
    # Determine organism name for folder naming
    local organism_name
    if [[ -n "$csv_organism_name" ]]; then
        # Use organism name extracted from CSV filename
        organism_name="$csv_organism_name"
    else
        # Fall back to organism name from database entry
        organism_name=$(execute_sqlite_cmd "$db_file" "SELECT organism FROM sra_entries WHERE sra_id = '$sra_id';" "true" 2>/dev/null)
        organism_name="${organism_name:-Unknown_Organism}"  # Default if not found
    fi
    local organism_safe=$(sanitize_filename "$organism_name")
    
    # Get species name from database for subfolder naming
    local species_name
    species_name=$(execute_sqlite_cmd "$db_file" "SELECT organism FROM sra_entries WHERE sra_id = '$sra_id';" "true" 2>/dev/null)
    species_name="${species_name:-Unknown_Species}"  # Default if not found
    local species_safe=$(sanitize_filename "$species_name")
    
    # Create new folder structure: organism/species/reads and organism/species/alignments
    local organism_dir="${organism_safe}"
    local species_dir="${organism_dir}/${species_safe}"
    local reads_dir="${species_dir}/reads"
    local alignments_dir="${species_dir}/alignments"
    local temp_dir="./tmp_processing/${sra_id}"

    # Cleanup any existing cache files before starting
    cleanup_sra_cache "$sra_id"

    # reset status to pending unless finished
    execute_sqlite_cmd "$db_file" "UPDATE sra_entries SET status = 'pending' WHERE status != 'finished';" >/dev/null
    # Create only temporary processing directory up front
    # NOTE: Organism/species output directories are created lazily
    #       ONLY when a sample passes the coverage threshold
    mkdir -p "$temp_dir"

    # Get technology type
    local technology
    debug_log "[Slot $slot] Retrieving technology type for $sra_id..."
    technology=$(execute_sqlite_cmd "$db_file" "SELECT technology FROM sra_entries WHERE sra_id = '$sra_id';")
    debug_log "[Slot $slot] Retrieved technology type: '$technology'"

    if [[ -z "$technology" ]]; then
        error_log "[Slot $slot] Failed to retrieve technology type for $sra_id"
        execute_sqlite_cmd "$db_file" "UPDATE sra_entries SET status = 'failed', message = 'Failed to retrieve technology type' WHERE sra_id = '$sra_id';" >/dev/null
        cleanup_sra_cache "$sra_id"
        rm -rf "$temp_dir"
        return 1
    fi

    # Set initial status to processing
    execute_sqlite_cmd "$db_file" "UPDATE sra_entries SET status = 'processing' WHERE sra_id = '$sra_id';" >/dev/null
    debug_log "[Slot $slot] Starting processing for $sra_id..."

    # Step 1: Prefetch
    prefetch --max-size 100G "$sra_id" >/dev/null 2>&1 &
    local prefetch_pid=$!
    track_subprocess "$prefetch_pid" "prefetch_${sra_id}"
    wait "$prefetch_pid"
    local prefetch_status=$?
    untrack_subprocess "$prefetch_pid"

    if [[ $prefetch_status -ne 0 ]]; then
        error_log "[Slot $slot] Failed to download $sra_id using prefetch"
        execute_sqlite_cmd "$db_file" "UPDATE sra_entries SET status = 'failed', message = 'Prefetch failed' WHERE sra_id = '$sra_id';" >/dev/null
        cleanup_sra_cache "$sra_id"
        rm -rf "$temp_dir"
        return 1
    fi

    # Step 2: Extract
    debug_log "[Slot $slot] Starting extraction for $sra_id..."
    local lock_file="/tmp/sra_${sra_id}.lock"
    (
        if ! flock -n 200; then
            debug_log "Waiting for lock on $sra_id..."
            flock 200
        fi

        fasterq-dump "$sra_id" \
            -O "$temp_dir" \
            -t "$temp_dir" \
            --split-files \
            -e "$threads" \
            >/dev/null 2>"${temp_dir}/error.log" &
        local dump_pid=$!
        track_subprocess "$dump_pid" "fasterq-dump_${sra_id}"
        wait "$dump_pid"
        local dump_status=$?
        untrack_subprocess "$dump_pid"
        exit "$dump_status"
    ) 200>"$lock_file"
    local extract_status=$?
    rm -f "$lock_file"

    if [[ $extract_status -ne 0 ]]; then
        error_log "[Slot $slot] Failed to extract $sra_id"
        execute_sqlite_cmd "$db_file" "UPDATE sra_entries SET status = 'failed', message = 'Extraction failed' WHERE sra_id = '$sra_id';" >/dev/null
        cleanup_sra_cache "$sra_id"
        rm -rf "$temp_dir"
        return 1
    fi

    # Step 3: Map to reference
    debug_log "[Slot $slot] Starting mapping for $sra_id..."
    
    # Find input files
    local input_files=("$temp_dir/${sra_id}"*.fastq)
    
    # Check if files exist
    if [[ ${#input_files[@]} -eq 0 ]] || [[ ! -f "${input_files[0]}" ]]; then
        error_log "[Slot $slot] No input files found for $sra_id"
        execute_sqlite_cmd "$db_file" "UPDATE sra_entries SET status = 'failed', message = 'No FASTQ files found' WHERE sra_id = '$sra_id';" >/dev/null
        cleanup_sra_cache "$sra_id"
        rm -rf "$temp_dir"
        return 1
    fi

    # Create coverage pipe
    local coverage_pipe="${temp_dir}/coverage_pipe"
    mkfifo "$coverage_pipe"

    # Start coverage calculation
    (
        local total_bases=0
        while read -r chrom pos depth; do
            ((total_bases += depth))
        done < "$coverage_pipe"

        # Calculate reference length
        local ref_length=$(samtools faidx "$reference" 2>/dev/null && awk '{sum+=$2} END {print sum}' "${reference}.fai")
        
        if [[ -n "$ref_length" ]] && [[ "$ref_length" -gt 0 ]]; then
            local coverage=$(echo "scale=2; $total_bases / $ref_length" | bc)
            coverage="${coverage:-0.00}"
            
            info_log "[Slot $slot] Coverage for $sra_id: ${coverage}x"
            
            if (( $(echo "$coverage > 1" | bc -l) )); then
                # Get organism name from database
                local species_name
                species_name=$(execute_sqlite_cmd "$db_file" "SELECT organism FROM sra_entries WHERE sra_id = '$sra_id';")
                species_name="${species_name:-Unknown Organism}"  # Default if not found
                pb push "VLE detected! $sra_id ($species_name) has ${coverage}x coverage"
                
                # Lazily create output directories ONLY after threshold is met
                mkdir -p "$organism_dir" "$species_dir" "$reads_dir" "$alignments_dir"
                
                # Compress and save reads
                for fastq_file in "$temp_dir"/*.fastq; do
                    local base_name=$(basename "$fastq_file")
                    pigz -c "$fastq_file" > "$reads_dir/${base_name}.gz"
                    
                done
                
                # Save the alignment file (BAM) for mapped reads
                if [[ -f "${temp_dir}/${sra_id}_mapped.bam" ]]; then
                    local bam_filename="${sra_id}_mapped.bam"
                    cp "${temp_dir}/${sra_id}_mapped.bam" "$alignments_dir/${bam_filename}"
                    info_log "[Slot $slot] Saved alignment file: $alignments_dir/${bam_filename}"
                fi
                
                execute_sqlite_cmd "$db_file" "UPDATE sra_entries SET status = 'finished', message = 'Coverage: ${coverage}x' WHERE sra_id = '$sra_id';" >/dev/null
                info_log "[Slot $slot] $sra_id processed successfully with ${coverage}x coverage"
            else
                execute_sqlite_cmd "$db_file" "UPDATE sra_entries SET status = 'finished', message = 'No coverage' WHERE sra_id = '$sra_id';" >/dev/null
                info_log "[Slot $slot] $sra_id has no coverage, discarding reads"
            fi
        fi
    ) &
    local coverage_pid=$!

    # Map based on technology
    local redirect_output="/dev/null"
    [[ "$DEBUG_MODE" == "true" ]] && redirect_output="/dev/stderr"

    case "$technology" in
        "illumina")
            if [[ ${#input_files[@]} -eq 1 ]]; then
                bwa-mem2 mem -t "$threads" "$reference" "${input_files[0]}" 2>$redirect_output | \
                samtools view -@ "$threads" -b -F 4 | \
                tee >(samtools sort -@ "$threads" -o "${temp_dir}/${sra_id}_mapped.bam" 2>$redirect_output) | \
                samtools sort -@ "$threads" 2>$redirect_output | \
                samtools depth -a - > "$coverage_pipe"
            else
                bwa-mem2 mem -t "$threads" "$reference" "${input_files[@]}" 2>$redirect_output | \
                samtools view -@ "$threads" -b -F 4 | \
                tee >(samtools sort -@ "$threads" -o "${temp_dir}/${sra_id}_mapped.bam" 2>$redirect_output) | \
                samtools sort -@ "$threads" 2>$redirect_output | \
                samtools depth -a - > "$coverage_pipe"
            fi
            ;;
        "nanopore"|"iontorrent"|"454"|"capillary")
            minimap2 -ax map-ont -t "$threads" "$reference" "${input_files[@]}" 2>$redirect_output | \
            samtools view -@ "$threads" -b -F 4 | \
            tee >(samtools sort -@ "$threads" -o "${temp_dir}/${sra_id}_mapped.bam" 2>$redirect_output) | \
            samtools sort -@ "$threads" 2>$redirect_output | \
            samtools depth -a - > "$coverage_pipe"
            ;;
        "pacbio")
            minimap2 -ax map-pb -t "$threads" "$reference" "${input_files[@]}" 2>$redirect_output | \
            samtools view -@ "$threads" -b -F 4 | \
            tee >(samtools sort -@ "$threads" -o "${temp_dir}/${sra_id}_mapped.bam" 2>$redirect_output) | \
            samtools sort -@ "$threads" 2>$redirect_output | \
            samtools depth -a - > "$coverage_pipe"
            ;;
        *)
            error_log "[Slot $slot] Unsupported technology: $technology"
            kill $coverage_pid 2>/dev/null
            cleanup_sra_cache "$sra_id"
            rm -rf "$temp_dir"
            return 1
            ;;
    esac

    wait $coverage_pid
    local mapping_status=$?

    if [[ $mapping_status -ne 0 ]]; then
        error_log "[Slot $slot] Mapping failed for $sra_id"
        execute_sqlite_cmd "$db_file" "UPDATE sra_entries SET status = 'failed', message = 'Mapping failed' WHERE sra_id = '$sra_id';" >/dev/null
        cleanup_sra_cache "$sra_id"
        rm -rf "$temp_dir"
        return 1
    fi

    # Increment the counter atomically
    local current_count
    {
        flock -x 200
        current_count=$(<"$counter_file")
        ((current_count++))
        echo "$current_count" > "$counter_file"
    } 200>"${counter_file}.lock"

    info_log "[Slot $slot] Completed processing $sra_id ($current_count out of $total_entries)"

    # Clean up before exiting
    debug_log "[Slot $slot] Cleaning up SRA cache for $sra_id..."
    cleanup_sra_cache "$sra_id"
    rm -rf "$temp_dir"
    return 0
}

# Function to start downloads
start_downloads() {
    local db_file="$1"
    local max_concurrent="$2"
    local reference="$3"  # Optional reference file
    local csv_organism_name="$4"  # Optional organism name from CSV filename
    
    debug_log "Starting processing with $max_concurrent concurrent processes..."
    
    # Setup trap to cleanup on exit
    trap 'cleanup_on_exit "$db_file"' INT TERM EXIT
    
    # Get all pending entries, excluding SOLiD and other unsupported technologies
    local pending_entries
    pending_entries=$(sqlite3 "$db_file" "SELECT sra_id FROM sra_entries WHERE status = 'pending' AND technology NOT IN ('other', 'solid');")
    
    # Check if we have any pending entries
    if [[ -z "$pending_entries" ]]; then
        info_log "No pending entries to process."
        return 0
    fi
    
    # Convert to array
    readarray -t sra_ids <<< "$pending_entries"
    local total_entries=${#sra_ids[@]}
    
    info_log "Found $total_entries entries to process"
    
    # Track child processes
    declare -a pids=()
    
    # Create a counter file for tracking progress
    local counter_file=$(mktemp)
    echo "0" > "$counter_file"
    
    # Start initial batch of processes
    local active=0
    local next_index=0
    
    # Start initial batch of processes
    for ((i=1; i<=max_concurrent && i<=total_entries; i++)); do
        process_sra_entry "${sra_ids[$next_index]}" "$db_file" "$i" "$counter_file" "$total_entries" "$reference" "$csv_organism_name" &
        pid=$!
        pids+=("$pid")
        info_log "[Slot $i] Started processing ${sra_ids[$next_index]}   (PID: $pid)"
        ((active++))
        ((next_index++))
    done
    
    # Process remaining entries as slots become available
    while ((next_index < total_entries)); do
        # Check for completed processes
        for i in "${!pids[@]}"; do
            if ! kill -0 "${pids[$i]}" 2>/dev/null; then
                # Process completed, remove from array
                unset "pids[$i]"
                ((active--))
                
                # Start a new process if there are more entries
                if ((next_index < total_entries)); then
                    process_sra_entry "${sra_ids[$next_index]}" "$db_file" "$((i+1))" "$counter_file" "$total_entries" "$reference" "$csv_organism_name" &
                    pid=$!
                    pids+=("$pid")
                    info_log "[Slot $((i+1))] Started processing ${sra_ids[$next_index]}   (PID: $pid)"
                    ((active++))
                    ((next_index++))
                fi
            fi
        done
        
        # Compact the array
        pids=("${pids[@]}")
        
        # Brief pause to avoid excessive CPU usage
        sleep 2
    done
    
    info_log "All processes started. Waiting for completion..."
    
    # Wait for all remaining processes to complete
    for pid in "${pids[@]}"; do
        if kill -0 "$pid" 2>/dev/null; then
            wait "$pid"
        fi
    done
    
    # Clean up counter file
    rm -f "$counter_file"
    
    # Remove the trap
    trap - INT TERM EXIT
    
    # Get final counts from database
    local final_stats
    final_stats=$(execute_sqlite_cmd "$db_file" "
        SELECT 
            COUNT(CASE WHEN status = 'finished' THEN 1 END) as processed,
            COUNT(CASE WHEN status = 'failed' THEN 1 END) as failed
        FROM sra_entries;")
    
    IFS='|' read -r processed failed <<< "$final_stats"
    
    info_log "Processing summary:"
    info_log "Total entries: $total_entries"
    info_log "Successfully processed: $processed"
    info_log "Failed: $failed"
}

# Function to count lines in a growing file
monitor_progress() {
    local file="$1"
    local total="$2"
    local last_count=0
    
    while true; do
        if [[ -f "$file" ]]; then
            local current_count=$(wc -l < "$file")
            if ((current_count > last_count)); then
                echo -ne "\rFetched $((current_count-1)) of $total entries ($(( (current_count-1) * 100 / total ))%)..."
                last_count=$current_count
            fi
        fi
        sleep 1
    done
}

# Function to print database statistics
print_db_stats() {
    local db_file="$1"
    info_log "Database statistics for $db_file:"
    
    # Get counts with retries - make sure the order in the SELECT matches the order in the read statement
    debug_log "Executing database statistics query..."
    local stats
    stats=$(execute_sqlite_cmd "$db_file" "
        SELECT 
            COUNT(*) as total,
            COUNT(CASE WHEN status = 'pending' THEN 1 END) as pending,
            COUNT(CASE WHEN status = 'finished' THEN 1 END) as finished,
            COUNT(CASE WHEN status = 'failed' THEN 1 END) as failed,
            COUNT(CASE WHEN technology = 'illumina' THEN 1 END) as illumina,
            COUNT(CASE WHEN technology = 'nanopore' THEN 1 END) as nanopore,
            COUNT(CASE WHEN technology = 'pacbio' THEN 1 END) as pacbio,
            COUNT(CASE WHEN technology = 'iontorrent' THEN 1 END) as iontorrent,
            COUNT(CASE WHEN technology = '454' THEN 1 END) as ls454,
            COUNT(CASE WHEN technology = 'capillary' THEN 1 END) as capillary,
            COUNT(CASE WHEN technology = 'solid' THEN 1 END) as solid,
            COUNT(CASE WHEN technology = 'other' THEN 1 END) as other
        FROM sra_entries;" "true" 2>/dev/null)
    
    debug_log "Raw stats output: $stats"
    
    # Run a simpler, separate query for each statistic to ensure accuracy
    local total=$(execute_sqlite_cmd "$db_file" "SELECT COUNT(*) FROM sra_entries;" "true" 2>/dev/null)
    local pending=$(execute_sqlite_cmd "$db_file" "SELECT COUNT(*) FROM sra_entries WHERE status = 'pending';" "true" 2>/dev/null)
    local finished=$(execute_sqlite_cmd "$db_file" "SELECT COUNT(*) FROM sra_entries WHERE status = 'finished';" "true" 2>/dev/null)
    local failed=$(execute_sqlite_cmd "$db_file" "SELECT COUNT(*) FROM sra_entries WHERE status = 'failed';" "true" 2>/dev/null)
    
    local illumina=$(execute_sqlite_cmd "$db_file" "SELECT COUNT(*) FROM sra_entries WHERE technology = 'illumina';" "true" 2>/dev/null)
    local nanopore=$(execute_sqlite_cmd "$db_file" "SELECT COUNT(*) FROM sra_entries WHERE technology = 'nanopore';" "true" 2>/dev/null)
    local pacbio=$(execute_sqlite_cmd "$db_file" "SELECT COUNT(*) FROM sra_entries WHERE technology = 'pacbio';" "true" 2>/dev/null)
    local iontorrent=$(execute_sqlite_cmd "$db_file" "SELECT COUNT(*) FROM sra_entries WHERE technology = 'iontorrent';" "true" 2>/dev/null)
    local ls454=$(execute_sqlite_cmd "$db_file" "SELECT COUNT(*) FROM sra_entries WHERE technology = '454';" "true" 2>/dev/null)
    local capillary=$(execute_sqlite_cmd "$db_file" "SELECT COUNT(*) FROM sra_entries WHERE technology = 'capillary';" "true" 2>/dev/null)
    local solid=$(execute_sqlite_cmd "$db_file" "SELECT COUNT(*) FROM sra_entries WHERE technology = 'solid';" "true" 2>/dev/null)
    local other=$(execute_sqlite_cmd "$db_file" "SELECT COUNT(*) FROM sra_entries WHERE technology = 'other';" "true" 2>/dev/null)

    # Display the statistics with default value of 0 if empty
    total="${total:-0}"
    pending="${pending:-0}"
    finished="${finished:-0}"
    failed="${failed:-0}"
    illumina="${illumina:-0}"
    nanopore="${nanopore:-0}"
    pacbio="${pacbio:-0}"
    iontorrent="${iontorrent:-0}"
    ls454="${ls454:-0}"
    capillary="${capillary:-0}"
    solid="${solid:-0}"
    other="${other:-0}"

    # Display the statistics
    info_log "  Total entries: $total"
    info_log "  Pending: $pending"
    info_log "  Finished: $finished"
    info_log "  Failed: $failed"
    info_log "  Illumina: $illumina"
    info_log "  Nanopore: $nanopore"
    info_log "  PacBio: $pacbio" 
    info_log "  IonTorrent: $iontorrent"
    info_log "  454: $ls454"
    info_log "  Capillary: $capillary"
    info_log "  SOLiD: $solid"
    info_log "  Other: $other"
}

# Function to clean up SRA cache files
cleanup_sra_cache() {
    local sra_id="$1"
    local cache_dir="/data/sra-cache/sra"
    
    debug_log "Cleaning up cache files for $sra_id..."
    
    # Clean up main SRA cache directory
    if [[ -d "$cache_dir" ]]; then
        rm -f "$cache_dir/${sra_id}.sra.lock" \
              "$cache_dir/${sra_id}.sra.prf" \
              "$cache_dir/${sra_id}.sra.tmp" \
              "$cache_dir/${sra_id}.sra" \
              "$cache_dir/${sra_id}.sra.vdbcache" 2>/dev/null
    fi
    
    # Clean up any temporary processing directories
    rm -rf "./tmp_processing/${sra_id}" 2>/dev/null
    
    debug_log "Cache cleanup completed for $sra_id"
    return 0
}

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to check required tools
check_required_tools() {
    local missing_tools=()
    
    # Check for SRA toolkit tools
    if ! command_exists prefetch || ! command_exists fasterq-dump; then
        missing_tools+=("SRA Toolkit (prefetch, fasterq-dump)")
    fi
    
    # Check for mapping tools
    if ! command_exists bwa-mem2; then
        missing_tools+=("bwa-mem2")
    fi

    if ! command_exists minimap2; then
        missing_tools+=("minimap2")
    fi
    
    # Check for samtools
    if ! command_exists samtools; then
        missing_tools+=("samtools")
    fi
    
    # Check for pigz
    if ! command_exists pigz; then
        missing_tools+=("pigz")
    fi
    
    # Report missing tools
    if [[ ${#missing_tools[@]} -gt 0 ]]; then
        error_log "Error: The following required tools are missing:"
        for tool in "${missing_tools[@]}"; do
            error_log "  - $tool"
        done
        error_log "Please install these tools before running this script with reference mapping."
        return 1
    fi
    
    return 0
}

# Function to generate BWA index
generate_bwa_index() {
    local reference="$1"
    
    # Check if index files already exist
    if [[ -f "${reference}.bwt" ]] && [[ -f "${reference}.pac" ]] && [[ -f "${reference}.ann" ]] && [[ -f "${reference}.amb" ]] && [[ -f "${reference}.sa" ]]; then
        debug_log "BWA index files already exist for $reference"
        return 0
    fi
    
    info_log "Generating BWA index for reference file: $reference"
    if ! bwa-mem2 index "$reference" 2>/dev/null; then
        error_log "Failed to generate BWA index for $reference"
        return 1
    fi
    
    debug_log "Successfully generated BWA index for $reference"
    return 0
}

# Function to check if an SRA entry exists in the database
check_sra_exists() {
    local db_file="$1"
    local sra_id="$2"
    
    # Check if the entry exists in the database
    local count
    count=$(sqlite3 "$db_file" "SELECT COUNT(*) FROM sra_entries WHERE sra_id = '$sra_id';")
    
    # Return true (0) if entry exists, false (1) otherwise
    return $((count == 0))
}

# Function to fetch data from SRA
fetch_sra_data() {
    local organism="$1"
    local csv_file="$2"
    local total_count="$3"
    local db_file="$4"
    local is_update="${5:-false}"  # New parameter to determine if this is an update operation
    local max_retries=5
    local retry_count=0
    local base_delay=10
    
    if [[ "$is_update" == "true" ]]; then
        # Use batch fetching for updates
        local batch_size=20
        local start=0
        local found_existing=false
        local temp_csv="${csv_file}.temp"
        local batch_file="${csv_file}.batch"
        
        # Write header first
        debug_log "Fetching header from SRA..."
        while ((retry_count < max_retries)); do
            if ((retry_count > 0)); then
                local delay=$((base_delay * 2 ** (retry_count - 1)))
                info_log "Retrying SRA query after ${delay}s delay (attempt $((retry_count+1))/$max_retries)..."
                sleep "$delay"
            fi
            
            if esearch -db sra -query "\"$organism\"[Organism] AND \"WGS\"[Strategy] AND \"GENOMIC\"[Source]" | \
               efetch -format runinfo -stop 1 > "$csv_file" 2>/dev/null; then
                break
            fi
            
            ((retry_count++))
            if ((retry_count == max_retries)); then
                error_log "Failed to fetch data from SRA after $max_retries attempts"
                return 1
            fi
        done
        
        # Reset retry count for batch fetching
        retry_count=0
        
        # Keep fetching batches until we find an existing entry or reach the end
        while ((start < total_count)) && [[ "$found_existing" == "false" ]]; do
            debug_log "Fetching batch starting at $start..."
            
            # Fetch current batch to temporary file with retries
            local batch_success=false
            retry_count=0
            
            while ((retry_count < max_retries)) && [[ "$batch_success" == "false" ]]; do
                if ((retry_count > 0)); then
                    local delay=$((base_delay * 2 ** (retry_count - 1)))
                    info_log "Retrying batch fetch after ${delay}s delay (attempt $((retry_count+1))/$max_retries)..."
                    sleep "$delay"
                fi
                
                if esearch -db sra -query "\"$organism\"[Organism] AND \"WGS\"[Strategy] AND \"GENOMIC\"[Source]" | \
                   efetch -format runinfo -start "$start" -stop "$((start + batch_size))" > "$batch_file" 2>/dev/null; then
                    batch_success=true
                    break
                fi
                
                ((retry_count++))
            done
            
            if [[ "$batch_success" == "false" ]]; then
                error_log "Failed to fetch batch after $max_retries attempts"
                break
            fi
            
            # Check if the batch file is empty or contains only a header
            if [[ ! -s "$batch_file" ]] || [[ $(wc -l < "$batch_file") -le 1 ]]; then
                debug_log "No more entries available after position $start"
                break
            fi
            
            # Process the batch file
            local found_in_batch=false
            while IFS= read -r line; do
                # Skip header line
                if [[ "$line" =~ ^Run, ]]; then
                    continue
                fi
                
                # Extract Run accession (first field)
                local sra_id=$(echo "$line" | cut -d',' -f1)
                
                if check_sra_exists "$db_file" "$sra_id"; then
                    debug_log "Found existing entry $sra_id, stopping batch fetch"
                    found_existing=true
                    found_in_batch=true
                    break
                else
                    # Append to main CSV if it's a new entry
                    echo "$line" >> "$csv_file"
                fi
            done < "$batch_file"
            
            # Break if we found an existing entry
            [[ "$found_in_batch" == "true" ]] && break
            
            ((start += batch_size))
            debug_log "Moving to next batch..."
            
            # Add a small delay between batches to avoid overwhelming the server
            sleep 2
        done
        
        # Clean up
        rm -f "$temp_csv" "$batch_file"
        
    else
        # Simple direct fetch for new database creation
        debug_log "Performing direct fetch from SRA..."
        while ((retry_count < max_retries)); do
            if ((retry_count > 0)); then
                local delay=$((base_delay * 2 ** (retry_count - 1)))
                info_log "Retrying SRA query after ${delay}s delay (attempt $((retry_count+1))/$max_retries)..."
                sleep "$delay"
            fi
            
            if esearch -db sra -query "\"$organism\"[Organism] AND \"WGS\"[Strategy] AND \"GENOMIC\"[Source]" | \
               efetch -format runinfo > "$csv_file" 2>/dev/null; then
                break
            fi
            
            ((retry_count++))
            if ((retry_count == max_retries)); then
                error_log "Failed to fetch data from SRA after $max_retries attempts"
                return 1
            fi
        done
    fi
    
    # Check if we got any entries
    local final_count=$(($(wc -l < "$csv_file")-1))
    if ((final_count > 0)); then
        info_log "Successfully fetched $final_count entries"
        return 0
    else
        info_log "No new entries were fetched"
        return 1
    fi
}

# Function to reset downloading status
reset_downloading_status() {
    local db_file="$1"
    info_log "Resetting downloading status..."
    execute_sqlite_cmd "$db_file" "UPDATE sra_entries SET status = 'pending' WHERE status = 'downloading';" >/dev/null
}

# Function to reset mapping status
reset_mapping_status() {
    local db_file="$1"
    info_log "Resetting mapping status..."
    execute_sqlite_cmd "$db_file" "UPDATE sra_entries SET status = 'pending' WHERE status = 'mapping';" >/dev/null
}

# Function to reset failed downloads
reset_failed_downloads() {
    local db_file="$1"
    local retry_all="$2"
    info_log "Resetting failed downloads..."
    if [[ "$retry_all" == "true" ]]; then
        execute_sqlite_cmd "$db_file" "UPDATE sra_entries SET status = 'pending', message = NULL WHERE status = 'failed';" >/dev/null
    else
        execute_sqlite_cmd "$db_file" "UPDATE sra_entries SET status = 'pending', message = NULL WHERE status = 'failed' AND message LIKE '%download%';" >/dev/null
    fi
}

# Function to cleanup on exit
cleanup_on_exit() {
    info_log "Cleaning up..."
    
    # Get all entries in processing state
    local db_file="$1"
    if [[ -n "$db_file" ]] && [[ -f "$db_file" ]]; then
        local processing_entries
        processing_entries=$(sqlite3 "$db_file" "SELECT sra_id FROM sra_entries WHERE status = 'processing';")
        
        if [[ -n "$processing_entries" ]]; then
            while IFS= read -r sra_id; do
                debug_log "Cleaning up SRA cache for $sra_id"
                cleanup_sra_cache "$sra_id"
            done <<< "$processing_entries"
            
            # Reset status to pending
            execute_sqlite_cmd "$db_file" "UPDATE sra_entries SET status = 'pending' WHERE status = 'processing';" >/dev/null
        fi
    fi
    
    # Cleanup subprocesses
    cleanup_subprocesses
    
    # Remove all temporary directories and files
    rm -rf "./tmp_processing"
    rm -f /tmp/sra_*.lock
    rm -f ./*.tmp
    rm -f ./*.lock
   
    # Note: {organism}/{species}/alignments directories are preserved intentionally
    # to keep the alignment files for successful mappings
    
    # Clean up SRA cache directory
    if [[ -d "/data/sra-cache/sra" ]]; then
        find "/data/sra-cache/sra" -type f -name "*.lock" -delete
        find "/data/sra-cache/sra" -type f -name "*.prf" -delete
        find "/data/sra-cache/sra" -type f -name "*.tmp" -delete
        find "/data/sra-cache/sra" -type f -name "*.cache" -delete
        find "/data/sra-cache/sra" -type f -name "*.vdbcache" -delete
    fi
    
    exit
}

# Main execution
main() {
    local csv_file=""
    local db_file=""
    local organism="fungi"
    local max_concurrent=5  # Default to 5 concurrent connections
    local retry_failed="false"
    local retry_all="false"
    local reference_file=""
    
    # Set up the main trap for cleanup
    trap 'cleanup_on_exit "$db_file"' INT TERM
    
    # Parse command line arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                print_usage
                exit 0
                ;;
            -c|--csv)
                if [[ -n $2 ]]; then
                    csv_file=$2
                    shift 2
                else
                    error_log "--csv requires a file path"
                    exit 1
                fi
                ;;
            -d|--db)
                if [[ -n $2 ]]; then
                    db_file=$2
                    shift 2
                else
                    error_log "--db requires a file path"
                    exit 1
                fi
                ;;
            -o|--organism)
                if [[ -n $2 ]]; then
                    organism=$2
                    shift 2
                else
                    error_log "--organism requires a value"
                    exit 1
                fi
                ;;
            -x|--connections)
                if [[ -n $2 ]]; then
                    max_concurrent=$2
                    shift 2
                fi
                ;;
            -r|--retry)
                retry_failed="true"
                shift
                ;;
            -R|--retry-all)
                retry_failed="true"
                retry_all="true"
                shift
                ;;
            -f|--reference)
                if [[ -n $2 ]]; then
                    reference_file=$2
                    shift 2
                else
                    error_log "--reference requires a file path"
                    exit 1
                fi
                ;;
            -D|--debug)
                DEBUG_MODE=true
                shift
                ;;
            *)
                error_log "Unknown argument: $1"
                print_usage
                exit 1
                ;;
        esac
    done
    
    # Check for reference file requirement
    if [[ -z "$reference_file" ]]; then
        error_log "ERROR: No reference file provided. A reference file is required for read mapping and coverage calculation."
        error_log "Please provide a reference file using the -f or --reference option."
        print_usage
        exit 1
    fi

    # Validate reference file
    if [[ ! -f "$reference_file" ]]; then
        error_log "ERROR: Reference file '$reference_file' does not exist"
        exit 1
    fi

    if [[ ! -r "$reference_file" ]]; then
        error_log "ERROR: Reference file '$reference_file' is not readable"
        exit 1
    fi

    # Check required tools and generate index for reference
    info_log "Checking for required tools..."
    if ! check_required_tools; then
        error_log "Missing required tools, cannot proceed with reference processing"
        exit 1
    fi

    info_log "Generating index for reference file..."
    if ! generate_bwa_index "$reference_file"; then
        error_log "Failed to generate index for reference file"
        exit 1
    fi
    
    [[ "$DEBUG_MODE" == "true" ]] && debug_log "Debug mode enabled"
    debug_log "Using $max_concurrent concurrent connections"
    
    # If database exists, process it
    if [[ -n "$db_file" ]]; then
        if [[ ! -f "$db_file" ]]; then
            error_log "Database file '$db_file' does not exist"
            exit 1
        fi
        if [[ ! -r "$db_file" ]]; then
            error_log "Database file '$db_file' is not readable"
            exit 1
        fi
        
        # Reset statuses if needed
        reset_downloading_status "$db_file"
        reset_mapping_status "$db_file"
        
        # Reset failed downloads if requested
        if [[ "$retry_failed" == "true" ]]; then
            reset_failed_downloads "$db_file" "$retry_all"
        fi
        
        # Print database statistics
        print_db_stats "$db_file"
        
        # If organism is specified, update database with new entries
        if [[ -n "$organism" ]] && [[ "$organism" != "fungi" ]]; then
            info_log "Fetching new data from SRA for organism: $organism..."
            local organism_safe=$(sanitize_filename "$organism")
            local temp_csv="${organism_safe}temp_sra_wgs.csv"
            
            # Get total count
            info_log "Querying SRA database..."
            local total_count
            total_count=$(esearch -db sra -query "\"$organism\"[Organism] AND \"WGS\"[Strategy] AND \"GENOMIC\"[Source]" | \
                         grep -o "<Count>[0-9]*</Count>" | sed 's/<Count>\([0-9]*\)<\/Count>/\1/')
            
            if [[ -z "$total_count" ]] || [[ "$total_count" -eq 0 ]]; then
                error_log "No entries found in SRA for organism: $organism"
                exit 1
            fi
            
            # Fetch new entries in batches
            if ! fetch_sra_data "$organism" "$temp_csv" "$total_count" "$db_file" "true"; then
                error_log "Failed to fetch new entries from SRA"
                rm -f "$temp_csv"
                exit 1
            fi
            
            # If we got any new entries, process them
            if [[ -f "$temp_csv" ]] && [[ $(wc -l < "$temp_csv") -gt 1 ]]; then
                info_log "Processing new entries..."
                if ! process_csv "$temp_csv" "$db_file"; then
                    error_log "Failed to process CSV file"
                    exit 1
                fi
                
                # Print updated statistics
                print_db_stats "$db_file"
                
                # Start concurrent processing for new entries
                start_downloads "$db_file" "$max_concurrent" "$reference_file" "$organism"
            else
                info_log "No new entries found"
                rm -f "$temp_csv"
                
                # Print current database statistics
                print_db_stats "$db_file"
                
                # Continue with existing entries
                info_log "Continuing with existing entries..."
                start_downloads "$db_file" "$max_concurrent" "$reference_file" "$organism"
            fi
        else
            # If no organism specified, just process existing entries
            # Try to extract organism name from database filename
            local db_basename=$(basename "$db_file")
            local db_organism_name=""
            if [[ "$db_basename" =~ ^(.+)_sra_wgs\.db$ ]]; then
                db_organism_name="${BASH_REMATCH[1]}"
                info_log "Extracted organism name from database filename: $db_organism_name"
            fi
            start_downloads "$db_file" "$max_concurrent" "$reference_file" "$db_organism_name"
        fi
        exit 0
    fi
    
    # Handle CSV file processing
    if [[ -n "$csv_file" ]]; then
        # Extract organism name from CSV filename if no organism specified
        local csv_organism_name=""
        if [[ -z "$organism" ]] || [[ "$organism" == "fungi" ]]; then
            # Extract organism name from CSV filename
            local csv_basename=$(basename "$csv_file")
            # Remove _sra_wgs.csv suffix if present, otherwise remove .csv
            if [[ "$csv_basename" =~ ^(.+)_sra_wgs\.csv$ ]]; then
                csv_organism_name="${BASH_REMATCH[1]}"
            elif [[ "$csv_basename" =~ ^(.+)\.csv$ ]]; then
                csv_organism_name="${BASH_REMATCH[1]}"
            else
                csv_organism_name="${csv_basename%.*}"  # fallback: remove any extension
            fi
            info_log "Extracted organism name from CSV filename: $csv_organism_name"
        else
            csv_organism_name="$organism"
        fi
        
        # Create database filename based on CSV or organism
        if [[ -z "$db_file" ]]; then
            # If CSV filename ends with _sra_wgs.csv, just replace with .db
            local csv_basename=$(basename "$csv_file")
            if [[ "$csv_basename" =~ ^(.+)_sra_wgs\.csv$ ]]; then
                db_file="${BASH_REMATCH[1]}_sra_wgs.db"
            else
                # Otherwise, use sanitized organism name with _sra_wgs.db
                local organism_safe=$(sanitize_filename "$csv_organism_name")
                db_file="${organism_safe}_sra_wgs.db"
            fi
        fi
        
        # Create and populate database from CSV
        init_database "$db_file"
        process_csv "$csv_file" "$db_file"
        
        # Print database statistics
        print_db_stats "$db_file"
        
        # Start concurrent processing
        start_downloads "$db_file" "$max_concurrent" "$reference_file" "$csv_organism_name"
        exit 0
    fi
    
    # Otherwise, fetch data for the specified organism
    local organism_safe=$(sanitize_filename "$organism")
    db_file="${organism_safe}_sra_wgs.db"
    csv_file="${organism_safe}_sra_wgs.csv"
    
    info_log "Fetching new data from SRA for organism: $organism..."
    
    # Get total count
    info_log "Querying SRA database..."
    local total_count
    total_count=$(esearch -db sra -query "\"$organism\"[Organism] AND \"WGS\"[Strategy] AND \"GENOMIC\"[Source]" | \
                 grep -o "<Count>[0-9]*</Count>" | sed 's/<Count>\([0-9]*\)<\/Count>/\1/')
    
    if [[ -z "$total_count" ]] || [[ "$total_count" -eq 0 ]]; then
        error_log "No entries found in SRA for organism: $organism"
        exit 1
    fi
    
    info_log "Found $total_count entries in SRA, starting download..."
    
    # Start progress monitoring in background
    monitor_progress "$csv_file" "$total_count" &
    monitor_pid=$!
    
    # Write header first
    esearch -db sra -query "\"$organism\"[Organism] AND \"WGS\"[Strategy] AND \"GENOMIC\"[Source]" | \
    efetch -format runinfo -stop 1 > "$csv_file"
    
    # Fetch data from SRA
    if ! esearch -db sra -query "\"$organism\"[Organism] AND \"WGS\"[Strategy] AND \"GENOMIC\"[Source]" | \
         efetch -format runinfo > "$csv_file"; then
        kill $monitor_pid 2>/dev/null
        error_log "Failed to fetch data from SRA"
        exit 1
    fi
    
    # Stop progress monitoring
    kill $monitor_pid 2>/dev/null
    wait $monitor_pid 2>/dev/null
    
    # Show final count
    local final_count=$(($(wc -l < "$csv_file")-1))
    info_log "Completed fetching $final_count of $total_count entries from SRA"
    
    # Create and populate database
    init_database "$db_file"
    process_csv "$csv_file" "$db_file"
    
    # Print database statistics
    print_db_stats "$db_file"
    
    # Start concurrent processing
    start_downloads "$db_file" "$max_concurrent" "$reference_file" "$organism"
}

# Run the script
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi