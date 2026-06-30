#!/bin/bash

# SRA Trawler ENA - European Nucleotide Archive Version
# Based on sra-trawler.sh but uses ENA APIs instead of NCBI
# Eliminates NCBI rate limiting and bulk metadata download issues

# Debug mode flag
DEBUG_MODE=false

# Persistent log file.
# Every info/warn/error/debug message is also appended here (with a timestamp
# and level) so a later session can read the file and see exactly what the
# script did, even after the terminal output is gone. Override the location by
# exporting LOG_FILE before running, e.g.  LOG_FILE=/path/to/run.log ./sra-trawler-ena.sh ...
LOG_FILE="${LOG_FILE:-sra-trawler-ena.log}"

# Append a single timestamped, plain-text (no color codes) line to LOG_FILE.
# Failures to write are ignored so logging never aborts the run.
_write_log() {
    local level="$1"; shift
    printf '%s [%-5s] [pid %d] %s\n' \
        "$(date '+%Y-%m-%d %H:%M:%S')" "$level" "$$" "$*" >> "$LOG_FILE" 2>/dev/null || true
}

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

# ENA API configuration
ENA_API_BASE="https://www.ebi.ac.uk/ena/portal/api"
ENA_RATE_LIMIT=45  # Stay under 50 req/sec limit
ENA_REQUEST_DELAY=$(echo "scale=3; 1/$ENA_RATE_LIMIT" | bc)

# Function to print debug messages
# Always recorded to LOG_FILE; only echoed to the terminal in debug mode.
debug_log() {
    _write_log "DEBUG" "$*"
    if [[ "$DEBUG_MODE" == "true" ]]; then
        echo -n "[DEBUG] $*" >&2
        echo >&2
    fi
}

# Function to print info messages
# NOTE: all log output goes to stderr, never stdout. Several functions return
# data to their caller via stdout command substitution (e.g. download_ena_fastq
# prints its list of downloaded file paths). If a log line went to stdout it
# would be captured into that return value and word-split into bogus tokens --
# this is exactly what made a recovered-from corrupted-gzip download warning
# pollute the file list and trigger a spurious "No FASTQ files found".
info_log() {
    _write_log "INFO" "$*"
    echo -e "${GREEN}[INFO]${NC}  $*" >&2
}

# Function to print warning messages
warn_log() {
    _write_log "WARN" "$*"
    echo -e "${YELLOW}[WARN]${NC}  $*" >&2
}

# Function to print error messages
error_log() {
    _write_log "ERROR" "$*"
    echo -e "${RED}[ERROR]${NC} $*" >&2
}

# Function to print usage information
print_usage() {
    cat << EOF
Usage: $(basename "$0") [options]

SRA Trawler ENA - European Nucleotide Archive Version
Uses ENA APIs for fast, reliable sequencing data discovery and download.

Options:
    -h, --help          Show this help message and exit
    -c, --csv FILE      Path to input CSV file with ENA metadata
    -d, --db FILE       Path to database file (optional)
                        If not provided, creates new database as {organism}_sra_wgs.db
    -o, --organism STR  Organism to search for (e.g., "fungi", "bacteria")
    -f, --reference FILE Path to reference genome FASTA file (required)
    -x, --connections INT Number of concurrent processing slots (default: 2)
                        ENA allows 50 req/sec, much higher than NCBI
    -N, --ncbi-slots INT Of the -x slots, how many should prefer NCBI prefetch
                        instead of ENA direct download (default: 0). Spreads
                        load across both providers; each entry falls back to the
                        other provider if its preferred one fails. NCBI is more
                        CPU/disk-heavy (prefetch + fasterq-dump), so fasterq-dump
                        threads are split across the NCBI slots. 0 = ENA-first,
                        SRA-fallback (original behaviour).
    -P, --pairs-only   When a run has a proper _1/_2 pair AND a bare singletons
                        file, map the pair only and drop the singletons. Default
                        is to map pairs PLUS singletons (single-end) and merge,
                        since the bare file can hold the bulk of the reads.
    -s, --source STR   ENA library_source to keep (default: GENOMIC). Use e.g.
                        TRANSCRIPTOMIC, METAGENOMIC, METATRANSCRIPTOMIC, or
                        ANY/ALL to disable the filter. Non-GENOMIC values are
                        tagged into the auto-generated DB name (e.g.
                        fungi_transcriptomic_sra_wgs.db) so sources never mix.
    -m, --min-coverage NUM Minimum coverage threshold (default: 1)
    -r, --retry        Retry transient/retryable failures (status 'failed':
                        network, server, disk, mapping). Permanent 'no_data'
                        entries are left untouched.
    -R, --retry-all    Retry ALL failures, including permanent 'no_data' entries
                        (reads missing/corrupt in both ENA and SRA)
    -n, --notify STR   Enable push notifications with custom message
    -D, --debug        Enable debug mode for verbose output

Examples:
    # Search ENA for fungi datasets
    $(basename "$0") -d ena_fungi.db -o "fungi" -f reference.fasta

    # Search for bacteria with notifications
    $(basename "$0") -d ena_bacteria.db -o "bacteria" -f reference.fasta -n "Bacteria"

    # Use existing CSV metadata
    $(basename "$0") -d ena.db -c ena_metadata.csv -f reference.fasta

    # Retry failed downloads (download-related only)
    $(basename "$0") -d ena.db -f reference.fasta -r

    # Retry ALL failed entries
    $(basename "$0") -d ena.db -f reference.fasta -R

    # Higher concurrency (ENA supports it!)
    $(basename "$0") -d ena_fungi.db -o "fungi" -f reference.fasta -x 20

Database Schema:
    Same as V1, but uses ENA accessions and metadata
    - ena_id: ENA run accession (primary key)
    - sample_name: Sample identifier
    - library_name: Library name
    - organism: Scientific organism name
    - platform: Sequencing platform
    - model: Instrument model
    - technology: Technology type (illumina, nanopore, etc.)
    - bases: Total bases in dataset
    - center_name: Sequencing center
    - message: Processing status/coverage info
    - status: Processing status:
        pending    - not yet processed (includes entries without an FTP link,
                     which are fetched via SRA prefetch)
        processing - currently being downloaded/mapped
        finished   - coverage computed (reads kept if >= min coverage)
        failed     - transient/retryable failure; reset with -r
        no_data    - permanent: reads missing/corrupt in ENA AND SRA; -R only
    - coverage: Numeric coverage value (REAL)
    - processed_date: Timestamp when the sample was processed (ISO format)
    - bioproject_id: ENA/NCBI BioProject accession (study_accession)
    - submission_date: Date the data was first made public (first_public)

ENA Advantages:
    ✅ 50 requests/second (vs NCBI's 3/sec)
    ✅ Direct organism queries (no 15GB bulk downloads)
    ✅ Reliable HTTP APIs (no curl 56 errors)
    ✅ Direct FASTQ downloads (no SRA Toolkit needed)
    ✅ Same data as NCBI (mirrored databases)

EOF
}

# Function to get taxonomy ID using NCBI E-utilities
lookup_taxonomy_id() {
    local organism="$1"

    debug_log "Looking up taxonomy ID for: $organism"

    # Try NCBI E-utilities taxonomy lookup
    local tax_response
    tax_response=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=${organism// /%20}" 2>/dev/null)

    if [[ -n "$tax_response" ]]; then
        # Extract first taxonomy ID from XML response
        local tax_id=$(echo "$tax_response" | grep -o '<Id>[0-9]\+</Id>' | head -1 | sed 's/<[^>]*>//g')
        if [[ -n "$tax_id" && "$tax_id" =~ ^[0-9]+$ ]]; then
            debug_log "Found taxonomy ID: $tax_id for $organism"
            echo "$tax_id"
            return 0
        fi
    fi

    debug_log "No taxonomy ID found for: $organism"
    echo ""
    return 1
}

# Function to get ENA taxonomy ID for organism (enhanced)
get_ena_taxonomy_id() {
    local organism="$1"

    # First try hardcoded common categories (fast)
    case "${organism,,}" in
        "fungi"|"fungus")
            echo "4751"  # Fungi superkingdom
            return 0
            ;;
        "bacteria"|"bacterial")
            echo "2"     # Bacteria superkingdom
            return 0
            ;;
        "archaea"|"archaeal")
            echo "2157"  # Archaea superkingdom
            return 0
            ;;
        "virus"|"viral")
            echo "10239" # Viruses superkingdom
            return 0
            ;;
        "plant"|"plants"|"plantae")
            echo "33090" # Viridiplantae
            return 0
            ;;
        "animal"|"animals"|"metazoa")
            echo "33208" # Metazoa
            return 0
            ;;
        # Common specific organisms
        "ascochyta rabiei")
            echo "5454"
            return 0
            ;;
        "saccharomyces cerevisiae")
            echo "4932"
            return 0
            ;;
        "escherichia coli")
            echo "511145"  # E. coli str. K-12 substr. MG1655
            return 0
            ;;
        "homo sapiens"|"human")
            echo "9606"
            return 0
            ;;
        *)
            # Try NCBI taxonomy lookup for specific species
            local tax_id=$(lookup_taxonomy_id "$organism")
            if [[ -n "$tax_id" ]]; then
                echo "$tax_id"
                return 0
            fi

            # Return empty to indicate organism name search should be used
            echo ""
            return 1
            ;;
    esac
}

# Function to sanitize filename
sanitize_filename() {
    local name="$1"
    # Convert to lowercase, replace spaces and special chars with underscore, then compress multiple underscores
    echo "$name" | tr '[:upper:]' '[:lower:]' | tr -c '[:alnum:]' '_' | tr -s '_' | sed 's/^_\|_$//g'
}

# Function to aggressively clean up temporary files for an entry
ena_cleanup() {
    local ena_id="$1"

    # Clean up temp processing directory immediately
    rm -rf "./tmp_processing/${ena_id}" 2>/dev/null

    # Clean up any temp files that might be left around in current directory
    rm -f "./tmp_ena_"* 2>/dev/null
    rm -f "./tmp_processing/${ena_id}"* 2>/dev/null

    # Clean up SRA cache for this entry
    cleanup_sra_cache "$ena_id"

    debug_log "Aggressive cleanup completed for $ena_id"
    return 0
}

# Function to clean up orphaned temp directories (older than 1 hour)
ena_cleanup_orphaned() {
    debug_log "Cleaning up orphaned temp directories and files..."

    # Remove temp processing directories older than 1 hour
    find "./tmp_processing" -maxdepth 1 -type d -name "ERR*" -o -name "SRR*" -o -name "DRR*" 2>/dev/null | \
    while read -r dir; do
        if [[ -d "$dir" ]] && [[ $(find "$dir" -maxdepth 0 -mmin +60 2>/dev/null) ]]; then
            debug_log "Removing orphaned directory: $dir"
            rm -rf "$dir" 2>/dev/null
        fi
    done

    # Remove orphaned temp files in current directory older than 1 hour
    find . -maxdepth 1 -name "tmp_ena_*" -type f -mmin +60 -delete 2>/dev/null

    return 0
}

# Function to clean up all temp files in current directory
ena_cleanup_all_temp() {
    debug_log "Cleaning up all temp files in current directory..."

    # Remove all temp files created by this script
    rm -f ./tmp_ena_* 2>/dev/null
    rm -f ./counter_ena_* 2>/dev/null
    rm -f ./sql_output_* 2>/dev/null
    rm -f ./sql_error_* 2>/dev/null

    # Remove tmp_processing directory
    rm -rf ./tmp_processing 2>/dev/null

    return 0
}

# Function to execute SQLite command with retries
execute_sqlite_cmd() {
    local db_file="$1"
    local cmd="$2"
    local max_retries=10
    local attempt=1
    local lock_file="${db_file}.lock"
    local result=""
    local suppress_debug="${3:-false}"

    # Ensure the database directory exists
    local db_dir=$(dirname "$db_file")
    if [[ ! -d "$db_dir" ]]; then
        mkdir -p "$db_dir" || {
            error_log "Failed to create database directory: $db_dir"
            return 1
        }
    fi

    while ((attempt <= max_retries)); do
        local tmp_output=$(mktemp ./sql_output_XXXX)
        local tmp_error=$(mktemp ./sql_error_XXXX)

        if [[ "$suppress_debug" != "true" ]]; then
            debug_log "Executing SQLite command (attempt $attempt/$max_retries)..."
        fi

        if (
            flock -x -w 10 200 || exit 1
            # Feed the SQL on stdin rather than as a command-line argument.
            # Large batches (e.g. 500 concatenated INSERTs) can exceed the OS
            # ARG_MAX limit, which made sqlite3 fail with
            # "Argument list too long" and silently drop the whole batch.
            sqlite3 "$db_file" > "$tmp_output" 2>"$tmp_error" <<<"$cmd"
        ) 200>"$lock_file"; then
            result=$(<"$tmp_output")
            rm -f "$tmp_output" "$tmp_error" "$lock_file"
            echo "$result"
            return 0
        else
            local error_msg=$(<"$tmp_error")

            if [[ "$suppress_debug" != "true" ]]; then
                debug_log "SQLite error (attempt $attempt): $error_msg"
            fi

            rm -f "$tmp_output" "$tmp_error"

            if ((attempt == max_retries)); then
                error_log "Failed to execute SQLite command after $max_retries attempts"
                error_log "Last error: $error_msg"
                rm -f "$lock_file"
                return 1
            fi

            # Random backoff between 1 and 10 seconds so that several
            # competing processes don't retry in lockstep
            local retry_delay=$((RANDOM % 10 + 1))
            if [[ "$suppress_debug" != "true" ]]; then
                debug_log "Retrying in ${retry_delay}s (attempt $attempt/$max_retries)..."
            fi
            sleep $retry_delay
        fi
        ((attempt++))
    done

    rm -f "$lock_file"
    return 1
}

# Function to initialize database schema
init_database() {
    local db_file="$1"

    info_log "Initializing ENA database schema..."

    if ! execute_sqlite_cmd "$db_file" "
        CREATE TABLE IF NOT EXISTS ena_entries (
            ena_id TEXT PRIMARY KEY,
            sample_name TEXT,
            library_name TEXT,
            organism TEXT,
            platform TEXT,
            model TEXT,
            technology TEXT,
            bases INTEGER,
            center_name TEXT,
            fastq_ftp TEXT,
            fastq_aspera TEXT,
            message TEXT,
            status TEXT DEFAULT 'pending',
            coverage REAL,
            processed_date TEXT,
            bioproject_id TEXT,
            submission_date TEXT
        );"; then
        error_log "Failed to initialize database schema"
        return 1
    fi

    if ! execute_sqlite_cmd "$db_file" "SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='ena_entries';" | grep -q "1"; then
        error_log "Failed to verify database initialization"
        return 1
    fi

    debug_log "Database initialized successfully"
    return 0
}

# Function to determine technology from platform and model
determine_technology() {
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

    # Ion Torrent platforms
    if [[ "$platform" == "ION_TORRENT" ]] ||
       [[ "$model" =~ (ION|PROTON|PGM|S5|GENEXUS) ]]; then
        echo "iontorrent"
        return 0
    fi

    # 454 pyrosequencing platforms
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

# Function to fetch ENA metadata for organism
# Function to backfill empty run-level center_name values from study-level data.
# ENA's read_run center_name is blank for many SRA-originated submissions, but
# the submitting institution is recorded on the parent study.
#
# We look up ONLY the specific studies referenced by empty-center rows, in small
# batched accession queries. This is deliberate: an unlimited taxon-wide study
# query (limit=0) gets truncated server-side -- especially right after the big
# metadata download -- returning a partial list (observed: 56, then 110 rows
# instead of ~41k) that silently backfills almost nothing. Bounded ~100-row
# responses don't get truncated, and querying by accession is also independent
# of which organism this run was invoked with. Idempotent and self-limiting:
# runs only when empty center_name rows exist (i.e. during a metadata fetch).
backfill_center_names_from_studies() {
    local db_file="$1"

    local empty_count
    empty_count=$(execute_sqlite_cmd "$db_file" "SELECT COUNT(*) FROM ena_entries WHERE center_name = '' OR center_name IS NULL;" "true" 2>/dev/null)
    if [[ -z "$empty_count" || "$empty_count" -eq 0 ]]; then
        debug_log "No empty center_name entries; skipping study-level backfill"
        return 0
    fi

    # Distinct studies that need a center_name.
    local bp_file
    bp_file=$(mktemp ./tmp_ena_bp_XXXX)
    execute_sqlite_cmd "$db_file" "SELECT DISTINCT bioproject_id FROM ena_entries WHERE (center_name = '' OR center_name IS NULL) AND bioproject_id IS NOT NULL AND bioproject_id != '';" "true" 2>/dev/null > "$bp_file"

    local bps
    mapfile -t bps < "$bp_file"
    rm -f "$bp_file"
    local study_total=${#bps[@]}
    if [[ "$study_total" -eq 0 ]]; then
        debug_log "No studies to look up for center_name backfill"
        return 0
    fi

    info_log "Backfilling center_name: looking up $study_total studies for $empty_count entries..."

    # Accumulate study_accession<TAB>center_name (non-empty only) for .import.
    local study_tsv
    study_tsv=$(mktemp ./tmp_ena_studies_XXXX)
    : > "$study_tsv"

    local batch_size=100 i j q resp
    for (( i=0; i<study_total; i+=batch_size )); do
        q=""
        for (( j=i; j<i+batch_size && j<study_total; j++ )); do
            [[ -z "${bps[j]}" ]] && continue
            [[ -n "$q" ]] && q="$q OR "
            q="${q}study_accession=\"${bps[j]}\""
        done
        [[ -z "$q" ]] && continue

        # POST (not GET) so the long OR-query isn't subject to URL-length limits.
        resp=$(curl -s -X POST "${ENA_API_BASE}/search" \
            --data-urlencode "query=$q" \
            --data-urlencode "result=study" \
            --data-urlencode "fields=study_accession,center_name" \
            --data-urlencode "format=tsv" 2>/dev/null)
        printf '%s\n' "$resp" | awk -F'\t' 'NR>1 && $2!="" {print}' >> "$study_tsv"

        (( (i / batch_size) % 10 == 0 )) && \
            debug_log "center_name backfill: queried $(( i+batch_size>study_total ? study_total : i+batch_size ))/$study_total studies"
    done

    local fetched_rows
    fetched_rows=$(wc -l < "$study_tsv")
    info_log "Retrieved center_name for $fetched_rows of $study_total studies"

    if [[ "$fetched_rows" -eq 0 ]]; then
        warn_log "No study-level center_name values retrieved; leaving center_name as-is"
        rm -f "$study_tsv"
        return 0
    fi

    # Load into a staging table and backfill in one sqlite3 session (.import and
    # the table must share a single connection). study_tsv has no header and
    # only non-empty centers, so we never overwrite an empty with another empty.
    sqlite3 "$db_file" <<EOF 2>/dev/null
DROP TABLE IF EXISTS _study_centers;
CREATE TABLE _study_centers (study_accession TEXT, center_name TEXT);
.mode tabs
.import '$study_tsv' _study_centers
CREATE INDEX _study_centers_idx ON _study_centers(study_accession);
UPDATE ena_entries
SET center_name = (
        SELECT s.center_name FROM _study_centers s
        WHERE s.study_accession = ena_entries.bioproject_id
        LIMIT 1)
WHERE (center_name = '' OR center_name IS NULL)
  AND bioproject_id IN (SELECT study_accession FROM _study_centers);
DROP TABLE _study_centers;
EOF

    local after
    after=$(execute_sqlite_cmd "$db_file" "SELECT COUNT(*) FROM ena_entries WHERE center_name = '' OR center_name IS NULL;" "true" 2>/dev/null)
    info_log "center_name backfill complete: filled $(( empty_count - after )), still empty ${after} (no study-level center)"

    rm -f "$study_tsv"
    return 0
}

fetch_ena_metadata() {
    local organism="$1"
    local db_file="$2"

    info_log "Fetching ENA metadata for organism: $organism"

    # Get taxonomy ID for the organism
    local tax_id=$(get_ena_taxonomy_id "$organism")

    # Build the optional library_source clause. Defaults to GENOMIC but is set by
    # -s/--source (SOURCE_FILTER); ANY/ALL disables the restriction. The clause is
    # appended to whichever query form we use below.
    local source_clause=""
    case "${SOURCE_FILTER^^}" in
        ""|"ANY"|"ALL")
            info_log "library_source filter: ANY (no source restriction)"
            ;;
        *)
            # URL-encode the value (spaces -> %20) and quote it if it contains
            # whitespace, e.g. library_source="VIRAL RNA".
            local src_enc="${SOURCE_FILTER// /%20}"
            [[ "$SOURCE_FILTER" == *" "* ]] && src_enc="%22${src_enc}%22"
            source_clause="%20AND%20library_source%3D${src_enc}"
            info_log "library_source filter: ${SOURCE_FILTER}"
            ;;
    esac

    # Build ENA query
    local query
    if [[ -n "$tax_id" ]]; then
        # Use taxonomy tree search (more accurate)
        query="tax_tree(${tax_id})%20AND%20library_strategy%3DWGS${source_clause}"
        info_log "Using taxonomy search for '$organism' (tax_id: $tax_id)"
    else
        # Fallback to organism name search with proper URL encoding
        warn_log "No taxonomy ID found for '$organism', using name search (less precise)"
        local encoded_organism=$(echo "$organism" | sed 's/ /%20/g')
        query="scientific_name%3D%22${encoded_organism}%22%20AND%20library_strategy%3DWGS${source_clause}"
        info_log "Using name search for '$organism'"
    fi

    # Fetch metadata from ENA Portal API using batching for large datasets
    local ena_url="${ENA_API_BASE}/search"
    local fields="run_accession,sample_accession,experiment_accession,study_accession,scientific_name,instrument_platform,instrument_model,library_name,library_strategy,library_source,library_selection,read_count,base_count,center_name,first_public,last_updated,experiment_title,study_title,study_alias,experiment_alias,fastq_bytes,fastq_md5,fastq_ftp,fastq_aspera,fastq_galaxy,submitted_bytes,submitted_md5,submitted_ftp,submitted_aspera,submitted_galaxy,sample_alias,sample_title,tax_id,sample_description"

    # ENA can return complete datasets in single requests (even 200k+ entries!)
    # Strategy: Try unlimited first, fall back to limited if that fails
    local temp_file=$(mktemp ./tmp_ena_metadata_XXXX)

    info_log "Fetching complete ENA metadata (attempting unlimited download)..."

    # First attempt: Get ALL data with no limit (limit=0 means unlimited)
    local params="query=${query}&result=read_run&fields=${fields}&format=tsv&limit=0"
    local full_url="${ena_url}?${params}"

    debug_log "Attempting unlimited request to get complete dataset"

    if curl -L -f -s "$full_url" > "$temp_file" 2>/dev/null; then
        local line_count=$(wc -l < "$temp_file")
        if [[ "$line_count" -ge 2 ]]; then
            info_log "Successfully downloaded complete dataset: $((line_count - 1)) entries"
        else
            error_log "Unlimited request returned no data"
            rm -f "$temp_file"
            return 1
        fi
    else
        # Unlimited request failed (maybe timeout/network), try limited fallback
        warn_log "Unlimited request failed, falling back to limited batch"
        local fallback_limit=100000

        params="query=${query}&result=read_run&fields=${fields}&format=tsv&limit=${fallback_limit}"
        full_url="${ena_url}?${params}"

        debug_log "Fallback request with limit $fallback_limit"

        if ! curl -L -f -s "$full_url" > "$temp_file"; then
            error_log "Fallback request also failed"
            rm -f "$temp_file"
            return 1
        fi

        local line_count=$(wc -l < "$temp_file")
        if [[ "$line_count" -lt 2 ]]; then
            error_log "Fallback request returned no data"
            rm -f "$temp_file"
            return 1
        fi

        if [[ "$line_count" -eq $(($fallback_limit + 1)) ]]; then
            warn_log "Hit fallback limit of $fallback_limit entries"
            warn_log "Complete dataset is larger - some entries may be missing"
            warn_log "Try running again or use more specific organism queries"
        fi

        info_log "Downloaded limited dataset: $((line_count - 1)) entries"
    fi

    # Check if we got any valid data
    local line_count=$(wc -l < "$temp_file")
    if [[ "$line_count" -lt 2 ]]; then
        error_log "No data returned from ENA for organism: $organism"
        rm -f "$temp_file"
        return 1
    fi

    info_log "Downloaded complete metadata: $((line_count - 1)) entries total"

    info_log "Processing $((line_count - 1)) ENA entries..."

    # Pre-filter to find only NEW entries (ULTRA FAST)
    info_log "Finding new entries to process..."
    local existing_ids_file=$(mktemp ./existing_ids_XXXX)
    local new_ids_file=$(mktemp ./new_ids_XXXX)

    # Get existing IDs (sorted)
    sqlite3 "$db_file" "SELECT ena_id FROM ena_entries;" | sort > "$existing_ids_file"
    local existing_count=$(wc -l < "$existing_ids_file")

    # Extract all IDs from TSV and find new ones
    tail -n +2 "$temp_file" | cut -f1 | sort > "$new_ids_file"
    local total_from_ena=$(wc -l < "$new_ids_file")

    # Find entries NOT in existing database
    local only_new_file=$(mktemp ./only_new_XXXX)
    comm -23 "$new_ids_file" "$existing_ids_file" > "$only_new_file"
    local new_count=$(wc -l < "$only_new_file")

    info_log "Existing entries: $existing_count, ENA total: $total_from_ena, New entries: $new_count"

    # Create filtered TSV with only new entries (MAJOR SPEEDUP)
    local filtered_tsv=$(mktemp ./filtered_tsv_XXXX)
    info_log "Pre-filtering TSV to only new entries..."

    # Extract header + lines for new entries only (much faster than processing all).
    # IMPORTANT: match the run accession on column 1 EXACTLY. The previous
    # `grep -F -f` did a substring match anywhere on the line, so a short
    # accession (e.g. ERR12345) would also pull in unrelated rows whose first
    # column merely contained it (e.g. ERR123456789), inflating counts and
    # muddying the "new entries" set. awk keys on $1 only.
    head -1 "$temp_file" > "$filtered_tsv"
    awk -F'\t' 'NR==FNR { ids[$1]=1; next } FNR>1 && ($1 in ids)' \
        "$only_new_file" "$temp_file" >> "$filtered_tsv"
    local filtered_count=$(( $(wc -l < "$filtered_tsv") - 1 ))

    rm -f "$new_ids_file" "$only_new_file"
    info_log "Filtered to $filtered_count new entries (skipping $((total_from_ena - filtered_count)) duplicates)"

    # Process the filtered TSV data (MUCH SMALLER NOW)
    local processed=0
    local skipped=0
    local current_batch=0
    local batch_size=500   # Balanced batch size (not too large to avoid SQL issues)
    local temp_sql=$(mktemp ./tmp_ena_sql_XXXX)

    # Technology counters
    local illumina_count=0
    local nanopore_count=0
    local pacbio_count=0
    local iontorrent_count=0
    local other_count=0

    # Count of rows whose fastq_ftp/fastq_aspera URLs arrived swapped and were
    # corrected by the host-based guard below.
    local swapped_urls=0

    # FTP-link counters: entries with an ENA FASTQ FTP link can be fetched
    # directly; entries WITHOUT one are still inserted as 'pending' and will be
    # retrieved later via the SRA prefetch fallback. Tracking these makes the
    # ENA-vs-SRA split visible at import time.
    local with_ftp_count=0
    local without_ftp_count=0

    # Start transaction
    echo "BEGIN TRANSACTION;" > "$temp_sql"

    # Skip header line and process data.
    # IMPORTANT: split on a NON-whitespace delimiter (SOH, \001), not the literal
    # tab. With `IFS=$'\t'`, tab is an IFS-whitespace character, so `read`
    # collapses runs of tabs and strips empty fields -- any blank column (e.g. an
    # empty library_name/instrument_model) shifts every later field left by one.
    # That mis-aligned library_strategy (so WGS rows looked non-WGS and were
    # skipped -> never inserted -> re-discovered as "new" on every run) and also
    # shifted fastq_ftp/fastq_aspera (the apparent URL "swap"). The input is
    # converted tab->\001 below; \001 never appears in ENA text, so empty fields
    # are preserved and columns stay aligned.
    while IFS=$'\001' read -r run sample experiment study scientific_name \
        instrument_platform instrument_model library_name library_strategy library_source \
        library_selection read_count base_count center_name first_public last_updated \
        experiment_title study_title study_alias experiment_alias fastq_bytes fastq_md5 \
        fastq_ftp fastq_aspera fastq_galaxy submitted_bytes submitted_md5 submitted_ftp \
        submitted_aspera submitted_galaxy sample_alias sample_title tax_id sample_description; do

        ((processed++))

        # Progress tracking (every 100 entries)
        if (( processed % 100 == 0 )); then
            info_log "Progress: $processed/$filtered_count entries processed ($(( processed * 100 / filtered_count ))%)"
        fi

        # All entries in filtered TSV are new - no duplicate checking needed!

        # Skip only entries without a run accession (cannot key the row).
        # Entries with an empty fastq_ftp are intentionally KEPT and inserted as
        # 'pending' so the processing loop can fetch them via the SRA fallback.
        # (Previously these were dropped without being inserted, so they were
        # re-discovered as "new" and re-skipped on every run.)
        if [[ -z "$run" ]]; then
            ((skipped++))
            continue
        fi

        # Skip if not WGS (the ENA query already filters to WGS, so this is
        # effectively never hit; kept as a defensive guard).
        if [[ "$library_strategy" != "WGS" ]]; then
            ((skipped++))
            continue
        fi

        # Skip rows whose library_source doesn't match the requested filter
        # (default GENOMIC; ANY/ALL disables this). The ENA query already filters
        # server-side, so this mainly guards older CSVs and rows with an empty
        # source field. For the GENOMIC default this drops TRANSCRIPTOMIC etc.,
        # whose coverage is expression-biased rather than genomic.
        case "${SOURCE_FILTER^^}" in
            ""|"ANY"|"ALL") ;;  # accept any source
            *)
                if [[ -n "$library_source" && "$library_source" != "${SOURCE_FILTER^^}" ]]; then
                    ((skipped++))
                    continue
                fi
                ;;
        esac

        # Defensive column-order guard for the FASTQ URLs. The downloader treats
        # fastq_ftp as an HTTP(S)-reachable URL (it prepends https:// and curls
        # it), so an Aspera "fasp:" URL landing there makes the download fail.
        # ~10% of the rows already in some DBs have fastq_ftp and fastq_aspera
        # swapped (cause unclear: older field order / an ENA hiccup). Route by
        # URL host instead of trusting the column position: the FTP/HTTP host
        # starts with "ftp", the Aspera host with "fasp". Only swap when it is
        # unambiguous (ftp field isn't ftp* but the aspera field is) so clean
        # rows are never disturbed.
        if [[ -n "$fastq_ftp" && "$fastq_ftp" != ftp* && "$fastq_aspera" == ftp* ]]; then
            local _swap_tmp="$fastq_ftp"
            fastq_ftp="$fastq_aspera"
            fastq_aspera="$_swap_tmp"
            ((swapped_urls++))
        fi

        # Inline sanitization (avoid 9 function calls per entry)
        local run_clean="${run//\'/\'\'}"
        local sample_clean="${sample//\'/\'\'}"
        local library_clean="${library_name//\'/\'\'}"
        local organism_clean="${scientific_name//\'/\'\'}"
        local platform_clean="${instrument_platform//\'/\'\'}"
        local model_clean="${instrument_model//\'/\'\'}"
        local center_clean="${center_name//\'/\'\'}"
        local fastq_ftp_clean="${fastq_ftp//\'/\'\'}"
        local fastq_aspera_clean="${fastq_aspera//\'/\'\'}"
        local study_clean="${study//\'/\'\'}"
        local first_public_clean="${first_public//\'/\'\'}"

        # Determine technology
        local tech=$(determine_technology "$instrument_platform" "$instrument_model")

        # Count technologies
        case "$tech" in
            "illumina") ((illumina_count++)) ;;
            "nanopore") ((nanopore_count++)) ;;
            "pacbio") ((pacbio_count++)) ;;
            "iontorrent") ((iontorrent_count++)) ;;
            *) ((other_count++)) ;;
        esac

        # Track whether this entry has a direct FASTQ FTP link. Either way it is
        # inserted as 'pending'; entries without a link fall through to the SRA
        # prefetch path during processing.
        if [[ -z "$fastq_ftp" ]]; then
            ((without_ftp_count++))
        else
            ((with_ftp_count++))
        fi

        # base_count must be inserted as a bare numeric literal (or NULL). ENA
        # occasionally returns an empty value; without this guard the empty
        # expansion produces invalid SQL (",,") that aborts the ENTIRE batch,
        # silently dropping up to batch_size new entries so they get
        # re-discovered as "new" on every subsequent run.
        local base_count_sql="$base_count"
        [[ "$base_count_sql" =~ ^[0-9]+$ ]] || base_count_sql="NULL"

        # Add to batch
        {
            echo "INSERT OR IGNORE INTO ena_entries (ena_id, sample_name, library_name, organism, platform, model, technology, bases, center_name, fastq_ftp, fastq_aspera, status, bioproject_id, submission_date)
                  VALUES ('$run_clean', '$sample_clean', '$library_clean', '$organism_clean', '$platform_clean', '$model_clean', '$tech', $base_count_sql, '$center_clean', '$fastq_ftp_clean', '$fastq_aspera_clean', 'pending', '$study_clean', '$first_public_clean');"
        } >> "$temp_sql"

        ((current_batch++))

        # Process batch if we've reached batch_size
        if ((current_batch >= batch_size)); then
            # Commit the current batch
            {
                echo "COMMIT;"
            } >> "$temp_sql"

            # Execute current batch with progress reporting.
            info_log "Executing SQL batch: $((processed - current_batch + 1))-$processed entries"
            local batch_output=$(execute_sqlite_cmd "$db_file" "$(cat $temp_sql)" "false" 2>&1)
            if [[ $? -ne 0 ]]; then
                # The batch could not be inserted; its entries are skipped this
                # run and will be re-discovered as "new" on the next fetch.
                error_log "Failed to process batch ending at entry $processed ($current_batch entries skipped)"
                error_log "SQL Error: $batch_output"
            else
                debug_log "Batch of $current_batch entries completed successfully (total: $processed)"
            fi

            # Start a fresh batch regardless of success/failure.
            current_batch=0
            echo "BEGIN TRANSACTION;" > "$temp_sql"
        fi
    done < <(tail -n +2 "$filtered_tsv" | tr '\t' '\001')

    # Process final batch if there are remaining entries
    if ((current_batch > 0)); then
        {
            echo "COMMIT;"
        } >> "$temp_sql"

        debug_log "Final batch has $current_batch entries"
        debug_log "Final batch SQL file size: $(wc -c < $temp_sql) bytes"
        debug_log "Final batch SQL (first 500 chars): $(head -c 500 $temp_sql)"

        # Try executing final batch with better error reporting
        info_log "Executing final SQL batch: $current_batch entries"
        local final_batch_output=$(sqlite3 "$db_file" < "$temp_sql" 2>&1)
        if [[ $? -ne 0 ]]; then
            error_log "Failed to process final batch"
            error_log "SQL Error: $final_batch_output"
            warn_log "Some entries in the final batch may not have been inserted"
            # Continue anyway - don't fail the entire operation
        else
            info_log "Final batch completed successfully"
        fi
    fi

    # Verify entries were inserted
    actual_count=$(execute_sqlite_cmd "$db_file" "SELECT COUNT(*) FROM ena_entries;" 2>/dev/null)

    if [[ -z "$actual_count" ]] || [[ "$actual_count" -eq 0 ]]; then
        error_log "Failed to insert entries into the database"
        rm -f "$temp_file" "$temp_sql"
        return 1
    fi

    # Cleanup
    rm -f "$temp_file" "$temp_sql" "$existing_ids_file" "$filtered_tsv"

    # Print summary
    info_log "Processing complete!"
    info_log "New entries processed: $processed"
    info_log "Entries skipped (no run accession / non-WGS): $skipped"
    info_log "New entries added to database: $((processed - skipped))"
    info_log "Download source split:"
    info_log "  With ENA FASTQ FTP link (direct download): $with_ftp_count"
    info_log "  Without FTP link (will use SRA prefetch):  $without_ftp_count"
    if (( swapped_urls > 0 )); then
        info_log "  Corrected swapped ftp/aspera URLs:         $swapped_urls"
    fi
    info_log "Technology distribution:"
    info_log "  Illumina: $illumina_count"
    info_log "  Nanopore: $nanopore_count"
    info_log "  PacBio: $pacbio_count"
    info_log "  IonTorrent: $iontorrent_count"
    info_log "  Other: $other_count"

    # ENA's run-level center_name is empty for ~10% of records even when the
    # submitting institution is recorded at the STUDY level. Backfill those by
    # looking up the specific studies the empty rows reference (batched, robust
    # to the truncation that a taxon-wide limit=0 query suffers).
    backfill_center_names_from_studies "$db_file"

    return 0
}

# Function to download FASTQ files from ENA
download_ena_fastq() {
    local ena_id="$1"
    local fastq_urls="$2"
    local output_dir="$3"

    debug_log "Downloading FASTQ files for $ena_id"

    # Split URLs by semicolon (ENA uses semicolon separator)
    local urls=(${fastq_urls//;/ })
    local downloaded_files=()

    for url in "${urls[@]}"; do
        if [[ -z "$url" ]]; then
            continue
        fi

        local filename=$(basename "$url")
        local output_path="$output_dir/$filename"

        debug_log "Downloading: $url"

        # Download with enhanced error reporting
        local max_retries=3
        local retry=0
        local success=false

        while [[ $retry -lt $max_retries ]] && [[ "$success" == "false" ]]; do
            # Make each (re)download attempt visible in the log so a recovered
            # transient failure (e.g. corrupted gzip) is traceable.
            debug_log "Download attempt $((retry + 1))/$max_retries for $filename"

            # Add protocol prefix if missing - use HTTPS for ENA URLs
            local full_url="$url"
            if [[ ! "$url" =~ ^https?:// ]] && [[ ! "$url" =~ ^ftp:// ]]; then
                full_url="https://$url"
            fi

            # Enhanced curl with error reporting
            local http_code
            local curl_exit_code

            # Get HTTP status and download file
            http_code=$(curl -L -f -s -w '%{http_code}' "$full_url" -o "$output_path")
            curl_exit_code=$?

            if [[ $curl_exit_code -eq 0 ]] && [[ -f "$output_path" ]]; then
                # Check file size and basic format validation
                local file_size=$(stat -f%z "$output_path" 2>/dev/null || stat -c%s "$output_path" 2>/dev/null || echo "0")

                if [[ $file_size -eq 0 ]]; then
                    ((retry++))
                    warn_log "Downloaded empty file for $filename (HTTP: $http_code, attempt $retry/$max_retries)"
                    echo "DOWNLOAD_ERROR:Empty file downloaded (HTTP: $http_code)" >&2
                    rm -f "$output_path"
                elif [[ ! "$filename" =~ \.fastq\.gz$ ]] && [[ ! "$filename" =~ \.fq\.gz$ ]]; then
                    # Check if it's actually a FASTQ file even if name doesn't match
                    local file_type=$(file -b "$output_path" 2>/dev/null || echo "unknown")
                    warn_log "Non-FASTQ filename detected: $filename (${file_size} bytes, type: $file_type)"
                    echo "DOWNLOAD_ERROR:Wrong file format - expected FASTQ.gz, got $file_type" >&2
                fi

                # Basic FASTQ format check for .gz files
                if [[ "$filename" =~ \.gz$ ]]; then
                    if ! gzip -t "$output_path" 2>/dev/null; then
                        ((retry++))
                        warn_log "Downloaded corrupted gzip file: $filename (HTTP: $http_code, attempt $retry/$max_retries)"
                        echo "DOWNLOAD_ERROR:Corrupted gzip file (HTTP: $http_code)" >&2
                        rm -f "$output_path"
                        continue
                    fi

                    # Check if it looks like FASTQ content (first few lines)
                    local first_char=$(zcat "$output_path" 2>/dev/null | head -c1)
                    if [[ "$first_char" != "@" ]] && [[ -n "$first_char" ]]; then
                        local actual_format=$(zcat "$output_path" 2>/dev/null | head -1 | cut -c1-10)
                        warn_log "File may not be FASTQ format: $filename starts with '$actual_format' (expected '@')"
                        echo "DOWNLOAD_ERROR:Not FASTQ format - starts with '$actual_format' (expected '@')" >&2
                    fi
                fi

                success=true
                downloaded_files+=("$output_path")
                debug_log "Downloaded: $filename (${file_size} bytes, HTTP: $http_code)"
            else
                ((retry++))
                # Detailed error reporting
                local error_detail="Unknown error"
                case $curl_exit_code in
                    6)  error_detail="Could not resolve host" ;;
                    7)  error_detail="Failed to connect" ;;
                    22) error_detail="HTTP error (status: $http_code)" ;;
                    23) error_detail="Write error (disk full?)" ;;
                    28) error_detail="Operation timeout" ;;
                    35) error_detail="SSL handshake failed" ;;
                    56) error_detail="Network receive error" ;;
                    *)  error_detail="Curl error $curl_exit_code (HTTP: $http_code)" ;;
                esac

                warn_log "Download failed for $filename: $error_detail (attempt $retry/$max_retries)"
                echo "DOWNLOAD_ERROR:$error_detail" >&2
                rm -f "$output_path"

                if [[ $retry -lt $max_retries ]]; then
                    sleep $((retry * 2))  # Exponential backoff
                fi
            fi
        done

        if [[ "$success" == "false" ]]; then
            # Store detailed error for database
            echo "DOWNLOAD_ERROR:Failed to download $filename after $max_retries attempts - check if file exists or format is correct" >&2
            error_log "Failed to download $filename after $max_retries attempts - check if file exists or format is correct"
            return 1
        fi

        # Rate limiting for ENA
        sleep "$ENA_REQUEST_DELAY"
    done

    # Return list of downloaded files
    printf '%s\n' "${downloaded_files[@]}"
    return 0
}

# Function to download and extract SRA files (fallback when no FASTQ URLs available)
download_sra_fastq() {
    local ena_id="$1"
    local output_dir="$2"
    # Thread budget for fasterq-dump. When several NCBI downloads run
    # concurrently (dual-source mode), the caller divides nproc among them to
    # avoid CPU oversubscription. Defaults to all cores when unset.
    local threads_arg="${3:-$(nproc)}"

    debug_log "Starting SRA download for $ena_id"

    # Step 1: Prefetch with size limit
    local prefetch_error_file="${output_dir}/${ena_id}_prefetch_error.log"
    local prefetch_output_file="${output_dir}/${ena_id}_prefetch_output.log"

    debug_log "Starting prefetch for $ena_id (max-size: 20G)"

    if prefetch --max-size 20G "$ena_id" >"$prefetch_output_file" 2>"$prefetch_error_file"; then
        debug_log "Prefetch completed for $ena_id"
    else
        local prefetch_error=$(cat "$prefetch_error_file" 2>/dev/null | head -1)
        error_log "Prefetch failed for $ena_id: $prefetch_error"
        echo "DOWNLOAD_ERROR:Prefetch failed - $prefetch_error" >&2
        rm -f "$prefetch_error_file" "$prefetch_output_file"
        return 1
    fi
    rm -f "$prefetch_error_file" "$prefetch_output_file"

    debug_log "Starting FASTQ extraction for $ena_id"

    # Step 2: Extract FASTQ files
    local dump_error_file="${output_dir}/${ena_id}_dump_error.log"
    local dump_output_file="${output_dir}/${ena_id}_dump_output.log"
    local threads="$threads_arg"

    # Use --split-3 (not --split-files): properly-paired reads go to _1/_2 and
    # reads missing their mate are written to a separate <acc>.fastq singletons
    # file. --split-files instead force-fills _1/_2, producing count-mismatched
    # pairs that break paired mapping/assembly. Long-read runs (nanopore/pacbio)
    # are unaffected -- they yield a single <acc>.fastq either way.
    if fasterq-dump "$ena_id" \
        -O "$output_dir" \
        -t "$output_dir" \
        --split-3 \
        -e "$threads" \
        >"$dump_output_file" 2>"$dump_error_file"; then

        debug_log "fasterq-dump completed for $ena_id"
    else
        local dump_output=$(cat "$dump_output_file" 2>/dev/null)
        local dump_error=$(cat "$dump_error_file" 2>/dev/null | head -1)
        error_log "Failed to extract FASTQ from $ena_id: $dump_error"
        debug_log "fasterq-dump failed stdout: $dump_output"
        echo "DOWNLOAD_ERROR:FASTQ extraction failed - $dump_error" >&2

        # Clean up SRA cache
        cleanup_sra_cache "$ena_id"
        rm -f "$dump_error_file" "$dump_output_file"
        return 1
    fi
    rm -f "$dump_error_file" "$dump_output_file"

    # Find and return the extracted FASTQ files
    local extracted_files=()

    # Find FASTQ files with flexible naming patterns
    local found_files=()

    debug_log "Looking for FASTQ files in: $output_dir"
    debug_log "Current working directory: $(pwd)"
    debug_log "Files in output_dir: $(ls -la "$output_dir" 2>/dev/null || echo 'directory not found')"

    # Look for any FASTQ files starting with the ena_id
    while IFS= read -r -d '' file; do
        found_files+=("$file")
    done < <(find "$output_dir" -maxdepth 1 -name "${ena_id}*.fastq" -print0 2>/dev/null)

    debug_log "Found files: ${found_files[*]}"

    # If no files found in output_dir, check current working directory (fasterq-dump sometimes ignores -O)
    if [[ ${#found_files[@]} -eq 0 ]]; then
        debug_log "No files in output_dir, checking current working directory..."
        while IFS= read -r -d '' file; do
            found_files+=("$file")
        done < <(find "." -maxdepth 1 -name "${ena_id}*.fastq" -print0 2>/dev/null)

        if [[ ${#found_files[@]} -gt 0 ]]; then
            debug_log "Found files in current directory: ${found_files[*]}"
            # Move files to output_dir
            for file in "${found_files[@]}"; do
                local basename=$(basename "$file")
                mv "$file" "$output_dir/$basename"
                debug_log "Moved $file to $output_dir/$basename"
            done
            # Update paths to reflect new location
            found_files=()
            while IFS= read -r -d '' file; do
                found_files+=("$file")
            done < <(find "$output_dir" -maxdepth 1 -name "${ena_id}*.fastq" -print0 2>/dev/null)
        fi
    fi

    if [[ ${#found_files[@]} -gt 0 ]]; then
        # Sort files to ensure consistent ordering (_1, _2, etc.)
        IFS=$'\n' extracted_files=($(sort <<<"${found_files[*]}"))
        debug_log "Found ${#extracted_files[@]} FASTQ files for $ena_id: ${extracted_files[*]##*/}"
    else
        error_log "No FASTQ files found after extraction for $ena_id in $output_dir"
        debug_log "Files in directory: $(ls -la "$output_dir" 2>/dev/null || echo 'directory not accessible')"
        echo "DOWNLOAD_ERROR:No FASTQ files found after SRA extraction" >&2
        return 1
    fi

    # Clean up SRA cache after successful extraction
    cleanup_sra_cache "$ena_id"

    debug_log "SRA download completed successfully for $ena_id (${#extracted_files[@]} files)"
    return 0
}

# Function to dynamically detect SRA cache directory
get_sra_cache_dir() {
    local sra_root=""
    local config_file="$HOME/.ncbi/user-settings.mkfg"

    debug_log "Detecting SRA cache directory..."

    # First, try to read from NCBI configuration file
    if [[ -f "$config_file" ]] && [[ -r "$config_file" ]]; then
        debug_log "Reading SRA configuration from: $config_file"
        # Extract repository root from config file
        sra_root=$(grep -E "^/repository/user/main/public/root\s*=" "$config_file" 2>/dev/null | \
                   sed -E 's|^/repository/user/main/public/root\s*=\s*"([^"]+)".*|\1|' | \
                   head -1)

        if [[ -n "$sra_root" ]] && [[ -d "$sra_root/sra" ]]; then
            debug_log "Found SRA cache from config: $sra_root/sra"
            echo "$sra_root/sra"
            return 0
        fi
    fi

    # Fallback locations to check
    local fallback_locations=(
        "/data/sra-cache/sra"
        "$HOME/.ncbi/public/sra"
        "/tmp/sra-cache/sra"
        "/var/cache/sra"
    )

    debug_log "Config file not found or invalid, trying fallback locations..."
    for location in "${fallback_locations[@]}"; do
        if [[ -d "$location" ]]; then
            debug_log "Found SRA cache at fallback location: $location"
            echo "$location"
            return 0
        fi
    done

    # If no existing directory found, use the most likely location
    debug_log "No existing SRA cache found, using default: /data/sra-cache/sra"
    echo "/data/sra-cache/sra"
    return 0
}

# Function to clean up SRA cache files
cleanup_sra_cache() {
    local sra_id="$1"
    local cache_dir="${SRA_CACHE_DIR:-$(get_sra_cache_dir)}"

    debug_log "Cleaning up SRA cache files for $sra_id in $cache_dir..."

    # Clean up main SRA cache directory
    if [[ -d "$cache_dir" ]]; then
        rm -f "$cache_dir/${sra_id}.sra.lock" \
              "$cache_dir/${sra_id}.sra.prf" \
              "$cache_dir/${sra_id}.sra.tmp" \
              "$cache_dir/${sra_id}.sra" \
              "$cache_dir/${sra_id}.sra.cache" \
              "$cache_dir/${sra_id}.sra.vdbcache" \
              2>/dev/null
    fi

    # Note: temp processing directories are cleaned up separately after file moves
    debug_log "SRA cache cleanup completed for $sra_id (temp processing dir preserved)"
    return 0
}

# Function to reset processing status entries
reset_processing_status() {
    local db_file="$1"
    info_log "Resetting processing status..."
    # Clear processed_date too: a 'pending' entry has not finished processing,
    # so it should not carry a (stale) processing timestamp.
    execute_sqlite_cmd "$db_file" "UPDATE ena_entries SET status = 'pending', processed_date = NULL WHERE status = 'processing';" >/dev/null
}

# Function to reset failed downloads
reset_failed_downloads() {
    local db_file="$1"
    local retry_all="$2"

    # Status semantics are now explicit, so no fragile message-keyword matching
    # is needed:
    #   'failed'  -> transient/retryable failure (network, server, disk, mapping)
    #   'no_data' -> permanent: reads missing or corrupt in both ENA and SRA
    if [[ "$retry_all" == "true" ]]; then
        info_log "Resetting ALL failed entries (transient + permanent) to pending..."
        local reset_count
        reset_count=$(execute_sqlite_cmd "$db_file" "UPDATE ena_entries SET status = 'pending', message = NULL, processed_date = NULL WHERE status IN ('failed', 'no_data'); SELECT changes();")
        info_log "Reset $reset_count entries to retry"
    else
        info_log "Resetting retryable (transient) failures to pending..."
        local reset_count
        reset_count=$(execute_sqlite_cmd "$db_file" "UPDATE ena_entries SET status = 'pending', message = NULL, processed_date = NULL WHERE status = 'failed'; SELECT changes();")
        info_log "Reset $reset_count transient failures to retry"
        info_log "(Permanent 'no_data' entries are left as-is; use -R to retry those too)"
    fi
}

# Comprehensive cleanup function for all scenarios
ena_full_cleanup() {
    debug_log "Performing full cleanup (temp files + SRA processes)..."

    # Kill any remaining SRA processes
    pkill -9 -f "prefetch" 2>/dev/null
    pkill -9 -f "fasterq-dump" 2>/dev/null

    # Clean up SRA cache directory
    local sra_cache_dir="${SRA_CACHE_DIR:-$(get_sra_cache_dir)}"
    if [[ -d "$sra_cache_dir" ]]; then
        find "$sra_cache_dir" -type f -name "*.lock" -delete 2>/dev/null
        find "$sra_cache_dir" -type f -name "*.prf" -delete 2>/dev/null
        find "$sra_cache_dir" -type f -name "*.tmp" -delete 2>/dev/null
        find "$sra_cache_dir" -type f -name "*.cache" -delete 2>/dev/null
        find "$sra_cache_dir" -type f -name "*.vdbcache" -delete 2>/dev/null
    fi

    # Clean up all temp files and directories
    ena_cleanup_all_temp
}

# Signal handlers for cleanup - will be set up later when actually processing

# Function to classify a download error message as either:
#   transient -> network/server hiccup that is worth retrying (-r resets these)
#   permanent -> the reads are genuinely missing or corrupt in BOTH ENA and the
#                SRA archive, so retrying will never help (only -R resets these)
#
# The download helpers emit "DOWNLOAD_ERROR:<detail>" lines; <detail> is what we
# classify here. When in doubt we err on the side of "transient" so a fixable
# glitch is never permanently dropped.
classify_download_error() {
    local msg="$1"

    # Permanent: accession/file does not exist, SRA yields no reads, or the
    # downloaded payload is not valid sequencing data. (No mid-pattern comments
    # here: a comment between '|'-joined case patterns is a syntax error.)
    case "$msg" in
        *"not found"* | *"cannot be found"* | *"does not exist"* | \
        *"no such"* | *"No such"* | *"404"* | \
        *"No FASTQ files found"* | *"FASTQ extraction failed"* | \
        *"Corrupted gzip"* | *"Not FASTQ format"* | *"Wrong file format"*)
            echo "permanent"
            return 0
            ;;
        *"Could not resolve host"* | *"Failed to connect"* | \
        *"timeout"* | *"Timeout"* | *"timed out"* | \
        *"SSL handshake"* | *"Network receive error"* | \
        *"Connection"* | *"connection"* | \
        *"Empty file"* | *"disk full"* | *"temporarily"* | \
        *" 500"* | *" 502"* | *" 503"* | *" 504"*)
            echo "transient"
            return 0
            ;;
        *)
            echo "transient"
            return 0
            ;;
    esac
}

# Function to process a single ENA entry
process_ena_entry() {
    local ena_id="$1"
    local db_file="$2"
    local reference="$3"
    local min_coverage="$4"
    local notify_string="$5"
    local slot="$6"
    local total_entries="$7"
    local counter_file="$8"
    local position="$9"
    local supplied_organism="${10}"   # organism from -o, used for the top-level folder
    local threads="$(nproc)"

    local progress="($position of $total_entries)"
    info_log "[Slot $slot] Processing $ena_id $progress..."

    # Fetch the per-sample organism (DB 'organism' column), technology and
    # download URLs in a single query.
    local metadata
    metadata=$(execute_sqlite_cmd "$db_file" "SELECT organism, technology, fastq_ftp, fastq_aspera FROM ena_entries WHERE ena_id = '$ena_id';" "true" 2>/dev/null)
    local species_name=$(echo "$metadata" | cut -d'|' -f1)
    local technology=$(echo "$metadata" | cut -d'|' -f2)
    local fastq_urls=$(echo "$metadata" | cut -d'|' -f3)
    local fastq_aspera_urls=$(echo "$metadata" | cut -d'|' -f4)
    species_name="${species_name:-Unknown_Species}"

    # Defensive guard for rows stored BEFORE the import-time fix, where the
    # fastq_ftp and fastq_aspera URLs are swapped. The download uses fastq_ftp
    # over HTTP(S), so a "fasp:" (Aspera) URL there would always fail. If the
    # ftp column isn't an ftp URL but the aspera column is, use the latter.
    if [[ -n "$fastq_urls" && "$fastq_urls" != ftp* && "$fastq_aspera_urls" == ftp* ]]; then
        debug_log "[Slot $slot] $ena_id has swapped ftp/aspera URLs; using ftp URL from aspera column"
        fastq_urls="$fastq_aspera_urls"
    fi

    # Folder layout: <supplied-organism>/<per-sample-organism>/{reads,alignments}
    #   top level = organism supplied via -o (the search query)
    #   sub level = this run's 'organism' column (per-sample scientific name)
    local organism_safe=$(sanitize_filename "${supplied_organism:-Unknown_Organism}")
    local species_safe=$(sanitize_filename "$species_name")
    local organism_dir="${organism_safe}"
    local species_dir="${organism_dir}/${species_safe}"
    local reads_dir="${species_dir}/reads"
    local alignments_dir="${species_dir}/alignments"
    local temp_dir="./tmp_processing/${ena_id}"

    # Create only temporary processing directory up front
    mkdir -p "$temp_dir"

    debug_log "[Slot $slot] Retrieved technology type: '$technology'"
    debug_log "[Slot $slot] FASTQ URLs: '$fastq_urls'"

    # Only fail if technology is missing - empty fastq_urls is OK for SRA fallback
    if [[ -z "$technology" ]]; then
        error_log "[Slot $slot] Failed to retrieve technology metadata for $ena_id"
        execute_sqlite_cmd "$db_file" "UPDATE ena_entries SET status = 'failed', message = 'Failed to retrieve technology type', processed_date = datetime('now') WHERE ena_id = '$ena_id';" >/dev/null
        increment_counter "$counter_file" "$ena_id" "$slot" "$total_entries"
        ena_cleanup "$ena_id"
        return 1
    fi

    # Log what we found
    if [[ -z "$fastq_urls" ]]; then
        debug_log "[Slot $slot] No FASTQ URLs found for $ena_id - will use SRA download"
    fi

    # Set initial status to processing
    execute_sqlite_cmd "$db_file" "UPDATE ena_entries SET status = 'processing' WHERE ena_id = '$ena_id';" >/dev/null
    debug_log "[Slot $slot] Starting processing for $ena_id..."

    # Download FASTQ files. Provider order is chosen per slot to spread load
    # across ENA (direct FASTQ) and NCBI (prefetch + fasterq-dump), then we fall
    # back to the OTHER provider if the first fails. Of the MAX_PARALLEL slots,
    # the highest-numbered NCBI_SLOTS prefer NCBI; the rest prefer ENA. With
    # NCBI_SLOTS=0 this is exactly the old behaviour (ENA first, SRA fallback).
    local downloaded_files=()
    local error_log_file="${temp_dir}/${ena_id}_download_error.log"
    local download_method=""
    local download_exit_code=1
    > "$error_log_file"

    # Per-slot provider preference.
    local primary secondary
    if (( NCBI_SLOTS > 0 )) && (( slot > MAX_PARALLEL - NCBI_SLOTS )); then
        primary="SRA";  secondary="ENA"
    else
        primary="ENA";  secondary="SRA"
    fi
    # Without an ENA FASTQ link the only option is NCBI (no ENA fallback).
    if [[ -z "$fastq_urls" ]]; then
        primary="SRA";  secondary=""
        debug_log "[Slot $slot] No FASTQ URLs for $ena_id - NCBI prefetch only"
    fi

    # Thread budget for any NCBI extraction: split cores across the NCBI slots
    # so concurrent fasterq-dump runs don't oversubscribe the CPU.
    local ncbi_threads="$threads"
    if (( NCBI_SLOTS > 0 )); then
        ncbi_threads=$(( threads / NCBI_SLOTS ))
        (( ncbi_threads < 1 )) && ncbi_threads=1
    fi

    local provider
    for provider in "$primary" "$secondary"; do
        [[ -z "$provider" ]] && continue
        (( download_exit_code == 0 )) && break   # already succeeded on primary

        if [[ "$provider" == "ENA" ]]; then
            [[ -z "$fastq_urls" ]] && continue
            debug_log "[Slot $slot] Trying ENA download for $ena_id"
            # download_ena_fastq returns the downloaded file paths on stdout;
            # its stderr (DOWNLOAD_ERROR:... detail) is kept for classification.
            downloaded_files=($(download_ena_fastq "$ena_id" "$fastq_urls" "$temp_dir" 2>>"$error_log_file"))
            if [[ ${#downloaded_files[@]} -gt 0 ]]; then
                download_exit_code=0
                download_method="ENA"
            else
                warn_log "[Slot $slot] ENA download failed for $ena_id"
                rm -f "${temp_dir}/"*.fastq* 2>/dev/null
                downloaded_files=()
            fi

        elif [[ "$provider" == "SRA" ]]; then
            debug_log "[Slot $slot] Trying NCBI prefetch download for $ena_id (threads: $ncbi_threads)"
            # Capture stderr (DOWNLOAD_ERROR:... detail) for accurate
            # classification, as the ENA path does.
            if download_sra_fastq "$ena_id" "$temp_dir" "$ncbi_threads" 2>>"$error_log_file"; then
                downloaded_files=()
                while IFS= read -r -d '' file; do
                    downloaded_files+=("$file")
                done < <(find "$temp_dir" -maxdepth 1 -name "${ena_id}*.fastq" -print0 2>/dev/null)

                if [[ ${#downloaded_files[@]} -gt 0 ]]; then
                    download_exit_code=0
                    download_method="SRA"
                else
                    debug_log "[Slot $slot] NCBI download succeeded but no FASTQ files found for $ena_id"
                    downloaded_files=()
                fi
            else
                warn_log "[Slot $slot] NCBI download failed for $ena_id"
                rm -f "${temp_dir}/"*.fastq* 2>/dev/null
            fi
        fi
    done

    # Check if any download method succeeded
    if [[ ${#downloaded_files[@]} -eq 0 ]] || [[ $download_exit_code -ne 0 ]]; then
        error_log "[Slot $slot] Both ENA and SRA download failed for $ena_id"

        # Extract the most recent detailed error captured from the ENA and/or
        # SRA download attempts. Both helpers emit "DOWNLOAD_ERROR:<detail>".
        local error_message="Download failed (no detail captured)"
        if [[ -f "$error_log_file" ]] && [[ -s "$error_log_file" ]]; then
            local detailed_error=$(grep -o "DOWNLOAD_ERROR:.*" "$error_log_file" | tail -1 | sed 's/DOWNLOAD_ERROR://')
            if [[ -n "$detailed_error" ]]; then
                error_message="$detailed_error"
            fi
        fi

        # Classify the failure so the status reflects whether it is worth
        # retrying. This applies whether the entry was downloaded via the ENA
        # FTP link or the SRA prefetch fallback:
        #   transient -> status 'failed'  (network/server/disk; -r retries it)
        #   permanent -> status 'no_data' (reads missing/corrupt everywhere;
        #                                  only -R will reset it)
        local failure_type
        failure_type=$(classify_download_error "$error_message")

        local new_status stored_message
        if [[ "$failure_type" == "permanent" ]]; then
            new_status="no_data"
            stored_message="Unavailable: ${error_message} (no retrievable reads in ENA or SRA)"
        else
            new_status="failed"
            stored_message="Retryable download error: ${error_message}"
        fi
        local stored_message_clean="${stored_message//\'/\'\'}"

        execute_sqlite_cmd "$db_file" "UPDATE ena_entries SET status = '$new_status', message = '$stored_message_clean', processed_date = datetime('now') WHERE ena_id = '$ena_id';" >/dev/null
        rm -f "$error_log_file"
        increment_counter "$counter_file" "$ena_id" "$slot" "$total_entries"
        ena_cleanup "$ena_id"
        return 1
    fi

    # Clean up error log file on success
    rm -f "$error_log_file"

    debug_log "[Slot $slot] Downloaded ${#downloaded_files[@]} FASTQ files for $ena_id using $download_method"

    # Use files directly - mapping tools handle .gz files natively
    local fastq_files=("${downloaded_files[@]}")
    debug_log "[Slot $slot] Starting mapping for $ena_id using ${#fastq_files[@]} FASTQ files"

    # Check if files exist
    if [[ ${#fastq_files[@]} -eq 0 ]] || [[ ! -f "${fastq_files[0]}" ]]; then
        error_log "[Slot $slot] No FASTQ files found for $ena_id"
        execute_sqlite_cmd "$db_file" "UPDATE ena_entries SET status = 'failed', message = 'No FASTQ files found', processed_date = datetime('now') WHERE ena_id = '$ena_id';" >/dev/null
        increment_counter "$counter_file" "$ena_id" "$slot" "$total_entries"
        ena_cleanup "$ena_id"
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

        local ref_length
        ref_length=$(samtools faidx "$reference" 2>/dev/null && awk '{sum+=$2} END {print sum}' "${reference}.fai")

        if [[ -n "$ref_length" ]] && [[ "$ref_length" -gt 0 ]]; then
            local coverage=$(echo "scale=2; $total_bases / $ref_length" | bc)
            coverage="${coverage:-0.00}"

            info_log "[Slot $slot] Coverage for $ena_id: ${coverage}x"

            # Check if coverage meets threshold
            if (( $(echo "$coverage >= $min_coverage" | bc -l) )); then
                # species_name is already set above and inherited by this subshell

                # Send push notification if enabled
                if [[ -n "$notify_string" ]]; then
                    pb push "${notify_string} detected! $ena_id ($species_name) has ${coverage}x coverage"
                fi

                # Create output directories only when threshold is met
                mkdir -p "$organism_dir" "$species_dir" "$reads_dir" "$alignments_dir"

                # Handle files based on download method (ENA vs SRA)
                debug_log "[Slot $slot] Coverage ${coverage}x >= ${min_coverage}x threshold - saving files for $ena_id"

                if [[ "$download_method" == "ENA" ]]; then
                    # ENA files are already .fastq.gz - just move them
                    debug_log "[Slot $slot] Moving ENA FASTQ files (already compressed) for $ena_id"
                    for fastq_file in "${fastq_files[@]}"; do
                        local base_name=$(basename "$fastq_file")
                        mv "$fastq_file" "$reads_dir/${base_name}"
                    done
                elif [[ "$download_method" == "SRA" ]]; then
                    # SRA files are .fastq - need compression before moving
                    debug_log "[Slot $slot] Compressing SRA FASTQ files with $COMPRESSION_TOOL -9 for $ena_id"
                    for fastq_file in "${fastq_files[@]}"; do
                        local base_name=$(basename "$fastq_file")
                        local compressed_file="$reads_dir/${base_name}.gz"

                        if [[ "$COMPRESSION_TOOL" == "pigz" ]]; then
                            if pigz -9 -c "$fastq_file" > "$compressed_file"; then
                                debug_log "[Slot $slot] Compressed with pigz: $compressed_file"
                            else
                                warn_log "[Slot $slot] pigz failed, trying gzip fallback for $fastq_file"
                                if gzip -9 -c "$fastq_file" > "$compressed_file"; then
                                    debug_log "[Slot $slot] Compressed with gzip fallback: $compressed_file"
                                else
                                    error_log "[Slot $slot] Both pigz and gzip failed for $fastq_file"
                                    continue
                                fi
                            fi
                        else
                            # Use gzip directly
                            if gzip -9 -c "$fastq_file" > "$compressed_file"; then
                                debug_log "[Slot $slot] Compressed with gzip: $compressed_file"
                            else
                                error_log "[Slot $slot] gzip failed for $fastq_file"
                                continue
                            fi
                        fi
                    done
                fi

                # Move BAM alignment file
                if [[ -f "${temp_dir}/${ena_id}_mapped.bam" ]]; then
                    mv "${temp_dir}/${ena_id}_mapped.bam" "$alignments_dir/${ena_id}_mapped.bam"
                    debug_log "[Slot $slot] Moved alignment file: $alignments_dir/${ena_id}_mapped.bam"
                fi

                execute_sqlite_cmd "$db_file" "UPDATE ena_entries SET status = 'finished', message = 'Coverage: above threshold', coverage = ${coverage}, processed_date = datetime('now') WHERE ena_id = '$ena_id';" >/dev/null
                info_log "[Slot $slot] $ena_id processed successfully with ${coverage}x coverage ($species_name)"
            else
                execute_sqlite_cmd "$db_file" "UPDATE ena_entries SET status = 'finished', message = 'Coverage: below threshold', coverage = ${coverage}, processed_date = datetime('now') WHERE ena_id = '$ena_id';" >/dev/null
                info_log "[Slot $slot] $ena_id has no coverage, discarding reads"
            fi
        fi
    ) &
    local coverage_pid=$!

    # Map based on technology (copied from V1 sra-trawler.sh)
    # Note: samtools -F 2308 filters out:
    #   - Flag 4: unmapped reads
    #   - Flag 256: secondary alignments (prevents NNNNNN sequences)
    #   - Flag 2048: supplementary alignments
    # This ensures only primary alignments are kept
    local redirect_output="/dev/null"
    [[ "$DEBUG_MODE" == "true" ]] && redirect_output="/dev/stderr"

    case "$technology" in
        "illumina")
            # A single run can expose up to three FASTQ files:
            #   <acc>_1 / <acc>_2 -> reads SRA treats as a proper pair
            #   <acc>.fastq(.gz)  -> reads SRA treats as unpaired/technical
            # These are DISJOINT subsets and EITHER can hold the bulk of the data
            # (the bare file is often the larger of the two), so we map both and
            # merge instead of guessing. Classify by filename suffix, NEVER by
            # download position: handing all files to one bwa call mis-pairs them
            # because bwa only consumes the first two read arguments as a pair.
            local r1="" r2="" single=""
            local f bn
            for f in "${fastq_files[@]}"; do
                bn=$(basename "$f")
                if   [[ "$bn" =~ _1\.(fastq|fq)(\.gz)?$ ]]; then r1="$f"
                elif [[ "$bn" =~ _2\.(fastq|fq)(\.gz)?$ ]]; then r2="$f"
                else single="$f"
                fi
            done

            local pair_ok="false"
            [[ -n "$r1" && -n "$r2" ]] && pair_ok="true"

            # Decide what to map:
            #   - If a proper _1/_2 pair exists, always map the pair. By default we
            #     ALSO map any bare singletons file (reads that lost their mate) as
            #     single-end and merge it into the same BAM, because that file can
            #     hold the bulk of the run's reads (bwa-mem2 supports both paired
            #     and single-end input). Pass --pairs-only to drop the singletons
            #     and screen on clean paired reads only (the previous behaviour).
            #   - Otherwise (true single-end run, or only an orphaned half-pair),
            #     fall back to single-end mapping of whatever reads exist.
            local single_files=()
            if [[ "$pair_ok" == "true" ]]; then
                if [[ -n "$single" ]]; then
                    if [[ "${PAIRS_ONLY:-false}" == "true" ]]; then
                        warn_log "[Slot $slot] $ena_id has paired + singletons files; --pairs-only set, dropping singletons ($(basename "$single"))"
                    else
                        info_log "[Slot $slot] $ena_id has paired + singletons files; mapping pairs plus singletons ($(basename "$single"))"
                        single_files+=("$single")
                    fi
                fi
            else
                [[ -n "$single" ]] && single_files+=("$single")
                [[ -n "$r1" ]] && single_files+=("$r1")
                [[ -n "$r2" ]] && single_files+=("$r2")
            fi

            # One coordinate-sorted BAM per source.
            local part_bams=()
            if [[ "$pair_ok" == "true" ]]; then
                bwa-mem2 mem -t "$threads" "$reference" "$r1" "$r2" 2>$redirect_output | \
                samtools view -@ "$threads" -b -F 2308 | \
                samtools sort -@ "$threads" -o "${temp_dir}/part_paired.bam" 2>$redirect_output
                [[ -s "${temp_dir}/part_paired.bam" ]] && part_bams+=("${temp_dir}/part_paired.bam")
            fi
            local idx=0
            for f in "${single_files[@]}"; do
                bwa-mem2 mem -t "$threads" "$reference" "$f" 2>$redirect_output | \
                samtools view -@ "$threads" -b -F 2308 | \
                samtools sort -@ "$threads" -o "${temp_dir}/part_single_${idx}.bam" 2>$redirect_output
                [[ -s "${temp_dir}/part_single_${idx}.bam" ]] && part_bams+=("${temp_dir}/part_single_${idx}.bam")
                ((idx++))
            done

            # Combine into the single mapped BAM the coverage stage expects, then
            # stream depth into the coverage pipe. Always write to the pipe (even
            # on failure) so the coverage subshell never blocks forever.
            if [[ ${#part_bams[@]} -eq 1 ]]; then
                mv "${part_bams[0]}" "${temp_dir}/${ena_id}_mapped.bam"
                samtools depth -a "${temp_dir}/${ena_id}_mapped.bam" > "$coverage_pipe"
            elif [[ ${#part_bams[@]} -gt 1 ]]; then
                samtools merge -f -@ "$threads" "${temp_dir}/${ena_id}_mapped.bam" "${part_bams[@]}" 2>$redirect_output
                samtools depth -a "${temp_dir}/${ena_id}_mapped.bam" > "$coverage_pipe"
            else
                error_log "[Slot $slot] No alignments produced for $ena_id"
                : > "$coverage_pipe"
            fi
            ;;
        "nanopore"|"iontorrent"|"454"|"capillary")
            minimap2 -ax map-ont -t "$threads" "$reference" "${fastq_files[@]}" 2>$redirect_output | \
            samtools view -@ "$threads" -b -F 2308 | \
            tee >(samtools sort -@ "$threads" -o "${temp_dir}/${ena_id}_mapped.bam" 2>$redirect_output) | \
            samtools sort -@ "$threads" 2>$redirect_output | \
            samtools depth -a - > "$coverage_pipe"
            ;;
        "pacbio")
            minimap2 -ax map-pb -t "$threads" "$reference" "${fastq_files[@]}" 2>$redirect_output | \
            samtools view -@ "$threads" -b -F 2308 | \
            tee >(samtools sort -@ "$threads" -o "${temp_dir}/${ena_id}_mapped.bam" 2>$redirect_output) | \
            samtools sort -@ "$threads" 2>$redirect_output | \
            samtools depth -a - > "$coverage_pipe"
            ;;
        *)
            error_log "[Slot $slot] Unsupported technology: $technology"
            kill $coverage_pid 2>/dev/null
            # Unsupported technology will never succeed on retry, so mark it
            # terminal ('no_data', only reset by -R) and stamp processed_date.
            # Previously this branch left the entry stuck in 'processing', so it
            # was reset to 'pending' and re-attempted on every run.
            local tech_clean="${technology//\'/\'\'}"
            execute_sqlite_cmd "$db_file" "UPDATE ena_entries SET status = 'no_data', message = 'Unsupported sequencing technology: ${tech_clean}', processed_date = datetime('now') WHERE ena_id = '$ena_id';" >/dev/null
            increment_counter "$counter_file" "$ena_id" "$slot" "$total_entries"
            ena_cleanup "$ena_id"
            return 1
            ;;
    esac

    wait $coverage_pid

    local mapping_status=$?
    debug_log "[Slot $slot] Mapping completed for $ena_id, status: $mapping_status"

    if [[ $mapping_status -ne 0 ]]; then
        error_log "[Slot $slot] Mapping failed for $ena_id"
        execute_sqlite_cmd "$db_file" "UPDATE ena_entries SET status = 'failed', message = 'Mapping failed', processed_date = datetime('now') WHERE ena_id = '$ena_id';" >/dev/null
        increment_counter "$counter_file" "$ena_id" "$slot" "$total_entries"
        ena_cleanup "$ena_id"
        return 1
    fi

    # FASTQ files and BAM files are handled in the coverage workflow above
    # Cleanup happens only after files are moved or coverage is insufficient

    # Increment the counter atomically (like V1)
    increment_counter "$counter_file" "$ena_id" "$slot" "$total_entries"

    # Clean up temporary directory aggressively (but keep temp files until final cleanup)
    debug_log "[Slot $slot] Processing completed successfully - temp files preserved for final cleanup"
    return 0
}

# Function to increment counter atomically (copied from V1)
increment_counter() {
    local counter_file="$1"
    local ena_id="$2"
    local slot="$3"
    local total_entries="$4"

    local current_count
    {
        flock -x 200
        current_count=$(<"$counter_file")
        ((current_count++))
        echo "$current_count" > "$counter_file"
    } 200>"${counter_file}.lock"

    # Cleanup temp files for this entry
    ena_cleanup "$ena_id"

    # Show completion message based on final status
    local final_status=$(execute_sqlite_cmd "$db_file" "SELECT status FROM ena_entries WHERE ena_id = '$ena_id';" | head -1)
    local final_message=$(execute_sqlite_cmd "$db_file" "SELECT message FROM ena_entries WHERE ena_id = '$ena_id';" | head -1)

    if [[ "$final_status" == "failed" ]]; then
        error_log "$ena_id: Failed - $final_message ($current_count out of $total_entries)"
    fi
    debug_log "[Slot $slot] Completed processing $ena_id ($current_count out of $total_entries)"
}

# Function to generate the BWA-MEM2 index for the reference, once.
# Without this, every `bwa-mem2 mem` call exits immediately printing nothing,
# which surfaces downstream as `[main_samview] fail to read the header from "-"`
# and "No alignments produced" for EVERY entry (typical on a fresh checkout where
# the index was never built). bwa-mem2 writes .0123/.amb/.ann/.bwt.2bit.64/.pac
# (NOT the classic BWA .bwt/.sa), so the existence check below tests for those.
generate_bwa_index() {
    local reference="$1"

    if [[ -f "${reference}.bwt.2bit.64" ]] && [[ -f "${reference}.0123" ]] && \
       [[ -f "${reference}.amb" ]] && [[ -f "${reference}.ann" ]] && \
       [[ -f "${reference}.pac" ]]; then
        debug_log "BWA-MEM2 index already present for $reference"
        return 0
    fi

    info_log "Generating BWA-MEM2 index for reference: $reference"
    local index_log="/dev/null"
    [[ "$DEBUG_MODE" == "true" ]] && index_log="/dev/stderr"
    if ! bwa-mem2 index "$reference" 2>"$index_log"; then
        error_log "Failed to generate BWA-MEM2 index for $reference"
        return 1
    fi

    debug_log "Successfully generated BWA-MEM2 index for $reference"
    return 0
}

# Function to start downloads with proper slot management (adapted from V1)
start_ena_downloads() {
    local db_file="$1"
    local max_concurrent="$2"
    local reference="$3"
    local min_coverage="$4"
    local notify_string="$5"
    local supplied_organism="$6"   # organism from -o, for the top-level output folder

    # Get pending entries that aren't 'other' or 'solid' technology (including entries without FASTQ URLs for SRA download)
    local pending_query="SELECT ena_id FROM ena_entries WHERE status = 'pending' AND technology NOT IN ('other', 'solid');"
    local pending_entries
    pending_entries=$(sqlite3 "$db_file" "$pending_query")

    # Check if we have any pending entries
    if [[ -z "$pending_entries" ]]; then
        info_log "No pending entries to process."
        return 0
    fi

    # Convert to array
    readarray -t ena_ids <<< "$pending_entries"
    local total_entries=${#ena_ids[@]}

    info_log "Found $total_entries entries to process"

    # Setup progress tracking
    local counter_file=$(mktemp ./counter_ena_XXXX)
    echo "0" > "$counter_file"

    # STABLE slot identity: slot_pid[s] is the PID currently running in slot s
    # (1..max_concurrent). When a slot's process finishes we relaunch the next
    # entry IN THE SAME SLOT NUMBER, so a slot's provider preference (ENA vs
    # NCBI) is fixed for the whole run.
    #
    # The previous implementation derived the slot number from a compacted array
    # index, which re-indexed every iteration. Because ENA slots finish quickly
    # and NCBI slots slowly, the NCBI tasks were repeatedly shifted to low
    # indices, so each freed NCBI slot was refilled with a low (ENA) slot number
    # -- and after the initial batch drained, every slot reverted to ENA.
    declare -a slot_pid=()
    local next_index=0 s

    # Start the initial batch: one entry per slot, slot number = s.
    for (( s=1; s<=max_concurrent && next_index<total_entries; s++ )); do
        process_ena_entry "${ena_ids[$next_index]}" "$db_file" "$reference" "$min_coverage" "$notify_string" "$s" "$total_entries" "$counter_file" "$((next_index + 1))" "$supplied_organism" &
        slot_pid[$s]=$!
        debug_log "[Slot $s] Started processing ${ena_ids[$next_index]}   (PID: ${slot_pid[$s]})"
        ((next_index++))
    done

    # Drain: when a slot's process exits, relaunch the next entry in that same
    # slot; if no entries remain, retire the slot. Loop until all slots retired.
    while (( ${#slot_pid[@]} > 0 )); do
        for s in "${!slot_pid[@]}"; do
            if ! kill -0 "${slot_pid[$s]}" 2>/dev/null; then
                if (( next_index < total_entries )); then
                    process_ena_entry "${ena_ids[$next_index]}" "$db_file" "$reference" "$min_coverage" "$notify_string" "$s" "$total_entries" "$counter_file" "$((next_index + 1))" "$supplied_organism" &
                    slot_pid[$s]=$!
                    debug_log "[Slot $s] Started processing ${ena_ids[$next_index]}   (PID: ${slot_pid[$s]})"
                    ((next_index++))
                else
                    unset "slot_pid[$s]"   # no more work for this slot
                fi
            fi
        done
        (( ${#slot_pid[@]} > 0 )) && sleep 2
    done

    info_log "All processes finished."
    wait 2>/dev/null   # reap any stragglers

    rm -f "$counter_file"
    return 0
}

# Function to check required tools
check_required_tools() {
    local missing_tools=()

    if ! command -v curl >/dev/null 2>&1; then
        missing_tools+=("curl")
    fi

    if ! command -v sqlite3 >/dev/null 2>&1; then
        missing_tools+=("sqlite3")
    fi

    if ! command -v bwa-mem2 >/dev/null 2>&1; then
        missing_tools+=("bwa-mem2")
    fi

    if ! command -v minimap2 >/dev/null 2>&1; then
        missing_tools+=("minimap2")
    fi

    if ! command -v samtools >/dev/null 2>&1; then
        missing_tools+=("samtools")
    fi

    if ! command -v pigz >/dev/null 2>&1; then
        missing_tools+=("pigz")
    fi

    if ! command -v bc >/dev/null 2>&1; then
        missing_tools+=("bc")
    fi

    if ! command -v prefetch >/dev/null 2>&1; then
        missing_tools+=("prefetch (SRA Toolkit)")
    fi

    if ! command -v fasterq-dump >/dev/null 2>&1; then
        missing_tools+=("fasterq-dump (SRA Toolkit)")
    fi

    # Check for compression tools (pigz preferred, gzip fallback)
    local compression_tool=""
    if command -v pigz >/dev/null 2>&1; then
        compression_tool="pigz"
        info_log "Using pigz for compression (best performance)"
    elif command -v gzip >/dev/null 2>&1; then
        compression_tool="gzip"
        warn_log "Using gzip fallback (pigz not available - install for better performance)"
    else
        missing_tools+=("pigz or gzip (compression)")
    fi
    export COMPRESSION_TOOL="$compression_tool"

    if [[ ${#missing_tools[@]} -gt 0 ]]; then
        error_log "Missing required tools:"
        for tool in "${missing_tools[@]}"; do
            error_log "  - $tool"
        done
        error_log "Please install missing tools and try again."
        return 1
    fi

    return 0
}

# Main function
main() {
    # Write a clear session banner so multiple runs are separable in the log,
    # and a later session can tell when this invocation started and with what
    # arguments.
    _write_log "INFO" "================ session start ================"
    _write_log "INFO" "Invocation: $(basename "$0") $*"
    _write_log "INFO" "Working dir: $(pwd)"
    info_log "Logging to: $LOG_FILE"

    local csv_file=""
    local db_file=""
    local organism=""
    local reference_file=""
    local max_parallel=2  # ENA allows higher concurrency
    local min_coverage=1
    local notify_string=""
    local retry_failed="false"
    local retry_all="false"
    local ncbi_slots=0    # of the -x slots, how many prefer NCBI prefetch (0 = ENA-first as before)
    local pairs_only="false"  # if true, drop singletons when a proper _1/_2 pair exists
    local source_filter="GENOMIC"  # ENA library_source filter; ANY/ALL disables it

    # Parse command line arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                print_usage
                exit 0
                ;;
            -c|--csv)
                csv_file="$2"
                shift 2
                ;;
            -d|--db)
                db_file="$2"
                shift 2
                ;;
            -o|--organism)
                organism="$2"
                shift 2
                ;;
            -f|--reference)
                reference_file="$2"
                shift 2
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
            -x|--connections)
                max_parallel="$2"
                shift 2
                ;;
            -N|--ncbi-slots)
                ncbi_slots="$2"
                shift 2
                ;;
            -P|--pairs-only)
                pairs_only="true"
                shift
                ;;
            -s|--source)
                source_filter="${2^^}"
                shift 2
                ;;
            -m|--min-coverage)
                min_coverage="$2"
                shift 2
                ;;
            -n|--notify)
                notify_string="$2"
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
    if [[ -z "$reference_file" ]]; then
        error_log "Reference file is required (-f/--reference)"
        exit 1
    fi

    # Validate / clamp the NCBI slot split and publish both values as globals so
    # process_ena_entry (run as a background job per slot) can read them.
    if ! [[ "$ncbi_slots" =~ ^[0-9]+$ ]]; then
        error_log "--ncbi-slots must be a non-negative integer (got: $ncbi_slots)"
        exit 1
    fi
    if (( ncbi_slots > max_parallel )); then
        warn_log "--ncbi-slots ($ncbi_slots) > --connections ($max_parallel); clamping to $max_parallel"
        ncbi_slots=$max_parallel
    fi
    NCBI_SLOTS=$ncbi_slots
    MAX_PARALLEL=$max_parallel
    PAIRS_ONLY=$pairs_only
    SOURCE_FILTER=$source_filter

    # Only require organism or CSV if database doesn't exist or is empty
    if [[ -z "$organism" ]] && [[ -z "$csv_file" ]]; then
        if [[ -z "$db_file" ]] || [[ ! -f "$db_file" ]]; then
            error_log "For new database: organism (-o/--organism) or CSV file (-c/--csv) is required"
            exit 1
        else
            # Check if database has entries to process
            local entry_count=$(sqlite3 "$db_file" "SELECT COUNT(*) FROM ena_entries;" 2>/dev/null || echo "0")
            if [[ "$entry_count" -eq 0 ]]; then
                error_log "Database is empty: organism (-o/--organism) or CSV file (-c/--csv) is required"
                exit 1
            fi
            info_log "Processing existing database with $entry_count entries (no new metadata fetch)"
        fi
    fi

    # When auto-naming the DB, tag it with the source unless it's the default
    # GENOMIC, so genomic and non-genomic data never share a default database
    # (the schema doesn't store library_source, so they can't be separated later).
    local db_source_tag=""
    case "$source_filter" in
        GENOMIC) ;;  # default -> keep the historical {organism}_sra_wgs.db name
        ANY|ALL) db_source_tag="_any" ;;
        *) db_source_tag="_$(echo "$source_filter" | tr 'A-Z ' 'a-z_')" ;;
    esac

    # Auto-generate database filename if not provided
    if [[ -z "$db_file" ]]; then
        if [[ -n "$csv_file" ]]; then
            # If CSV filename ends with _sra_wgs.csv, just replace with .db
            local csv_basename=$(basename "$csv_file")
            if [[ "$csv_basename" =~ ^(.+)_sra_wgs\.csv$ ]]; then
                db_file="${BASH_REMATCH[1]}_sra_wgs.db"
            else
                # Otherwise, use sanitized organism name with _sra_wgs.db
                local organism_safe=$(sanitize_filename "${organism:-ena}")
                db_file="${organism_safe}_sra_wgs.db"
            fi
        else
            # Use organism name (+ source tag) for database filename
            local organism_safe=$(sanitize_filename "$organism")
            db_file="${organism_safe}${db_source_tag}_sra_wgs.db"
        fi
        info_log "Auto-generated database filename: $db_file"
    fi

    if [[ ! -f "$reference_file" ]]; then
        error_log "Reference file does not exist: $reference_file"
        exit 1
    fi

    # Check for required tools
    if ! check_required_tools; then
        exit 1
    fi

    # Set up signal handlers for cleanup (only when actually processing, not for help)
    trap 'ena_full_cleanup' EXIT
    trap 'ena_full_cleanup; exit 130' INT
    trap 'ena_full_cleanup; exit 143' TERM

    # Detect and cache SRA directory once at startup
    export SRA_CACHE_DIR=$(get_sra_cache_dir)
    info_log "SRA cache directory: $SRA_CACHE_DIR"

    # Check for push notification tool
    if [[ -n "$notify_string" ]] && ! command -v pb >/dev/null 2>&1; then
        error_log "Push notification requested but 'pb' tool not found"
        error_log "Install pushbullet-bash or remove -n option"
        exit 1
    fi

    # Convert reference to absolute path
    reference_file=$(realpath "$reference_file")

    # Build the BWA-MEM2 index once up front (idempotent). Mapping fails silently
    # for every entry if this is missing, so do it before launching any slots.
    if ! generate_bwa_index "$reference_file"; then
        error_log "Cannot proceed without a BWA-MEM2 index for the reference"
        exit 1
    fi

    info_log "SRA Trawler ENA - Starting analysis"
    info_log "Database: $db_file"
    info_log "Reference: $reference_file"
    info_log "Max parallel: $max_parallel"
    if (( NCBI_SLOTS > 0 )); then
        info_log "Provider split: $(( max_parallel - NCBI_SLOTS )) ENA-first slot(s), $NCBI_SLOTS NCBI-first slot(s)"
    else
        info_log "Provider split: ENA-first with NCBI fallback (all slots)"
    fi
    info_log "Min coverage: ${min_coverage}x"

    # Initialize database
    if ! init_database "$db_file"; then
        exit 1
    fi

    # Reset statuses if needed (similar to V1 script)
    reset_processing_status "$db_file"

    # Reset failed downloads if requested
    if [[ "$retry_failed" == "true" ]]; then
        reset_failed_downloads "$db_file" "$retry_all"
    fi

    # Fetch metadata if organism specified
    if [[ -n "$organism" ]]; then
        info_log "Fetching ENA metadata for: $organism"
        if ! fetch_ena_metadata "$organism" "$db_file"; then
            error_log "Failed to fetch metadata for $organism"
            exit 1
        fi
    fi

    # TODO: Add CSV processing support if needed
    if [[ -n "$csv_file" ]]; then
        warn_log "CSV file processing not yet implemented"
    fi

    # Clean up any orphaned temp directories and files from previous runs
    ena_cleanup_orphaned
    ena_cleanup_all_temp

    # Cleanup already configured at script start

    # Top-level output folder is named after the supplied organism (-o). When
    # resuming an existing DB without -o, fall back to the organism encoded in
    # the database filename ({organism}_sra_wgs.db) so output still lands in a
    # sensibly-named folder.
    local folder_organism="$organism"
    if [[ -z "$folder_organism" ]]; then
        folder_organism=$(basename "$db_file")
        folder_organism="${folder_organism%_sra_wgs.db}"
        folder_organism="${folder_organism%.db}"
    fi

    # Process entries
    info_log "Starting parallel processing with $max_parallel slots..."
    info_log "Output folder root: $(sanitize_filename "$folder_organism")/"
    if ! start_ena_downloads "$db_file" "$max_parallel" "$reference_file" "$min_coverage" "$notify_string" "$folder_organism"; then
        error_log "Processing failed"
        exit 1
    fi

    # Final summary
    local total=$(execute_sqlite_cmd "$db_file" "SELECT COUNT(*) FROM ena_entries;" 2>/dev/null)
    local finished=$(execute_sqlite_cmd "$db_file" "SELECT COUNT(*) FROM ena_entries WHERE status = 'finished';" 2>/dev/null)
    local pending=$(execute_sqlite_cmd "$db_file" "SELECT COUNT(*) FROM ena_entries WHERE status = 'pending';" 2>/dev/null)
    local failed=$(execute_sqlite_cmd "$db_file" "SELECT COUNT(*) FROM ena_entries WHERE status = 'failed';" 2>/dev/null)
    local no_data=$(execute_sqlite_cmd "$db_file" "SELECT COUNT(*) FROM ena_entries WHERE status = 'no_data';" 2>/dev/null)

    info_log "Processing complete!"
    info_log "Total entries:        $total"
    info_log "Finished:             $finished"
    info_log "Pending (unprocessed): $pending"
    info_log "Failed (retryable, -r): $failed"
    info_log "No data (permanent, -R): $no_data"

    return 0
}

# Run the script
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi