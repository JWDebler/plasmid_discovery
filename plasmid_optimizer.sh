#!/usr/bin/env bash
set -euo pipefail

# Default values
THREADS=1

usage() {
        echo "Usage: $0 [--reference <ref.fa>] [--plasmids <dir>] [--output <dir>] [--threads <N>]" >&2
        echo "  --reference, -r      Reference plasmid FASTA file" >&2
        echo "  --plasmids, -p       Directory containing plasmid FASTA files" >&2
        echo "  --output, -o         Output directory for optimized plasmids (default: optimized_plasmids)" >&2
        echo "  --threads, -t        Number of threads (default: 1)" >&2
        echo "  --help, -h           Show this help" >&2
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
        case $1 in
                --reference|-r)
                        REFERENCE_ARG="$2"
                        shift 2
                        ;;
                --plasmids|-p)
                        PLASMIDS_DIR_ARG="$2"
                        shift 2
                        ;;
                --output|-o)
                        OUTPUT_DIR_ARG="$2"
                        shift 2
                        ;;
                --threads|-t)
                        THREADS="$2"
                        shift 2
                        ;;
                --help|-h)
                        usage
                        exit 0
                        ;;
                *)
                        echo "Unknown option: $1" >&2
                        usage
                        exit 1
                        ;;
        esac
done

# Validate required arguments
if [ -z "$REFERENCE_ARG" ] || [ -z "$PLASMIDS_DIR_ARG" ]; then
        echo "Error: Both --reference and --plasmids are required" >&2
        usage
        exit 1
fi

# Resolve paths
REFERENCE="$(cd "$(dirname "$REFERENCE_ARG")" && pwd)/$(basename "$REFERENCE_ARG")"
PLASMIDS_DIR="$(cd "$PLASMIDS_DIR_ARG" && pwd)"
OUTPUT_DIR="${OUTPUT_DIR_ARG:-optimized_plasmids}"

# Check if reference exists and is valid
if [ ! -f "$REFERENCE" ]; then
        echo "Error: Reference file not found: $REFERENCE" >&2
        exit 1
fi

if [ ! -d "$PLASMIDS_DIR" ]; then
        echo "Error: Plasmids directory not found: $PLASMIDS_DIR" >&2
        exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "==================================================================" >&2
echo "PLASMID OPTIMIZATION PIPELINE" >&2
echo "==================================================================" >&2
echo "Reference plasmid: $REFERENCE" >&2
echo "Plasmids directory: $PLASMIDS_DIR" >&2
echo "Output directory: $OUTPUT_DIR" >&2
echo "Threads: $THREADS" >&2
echo "==================================================================" >&2

# Function to get plasmid length
get_plasmid_length() {
        local fasta_file="$1"
        if command -v seqkit >/dev/null 2>&1; then
                seqkit fx2tab -nl "$fasta_file" | cut -f2
        else
                awk 'BEGIN{seq_len=0} /^>/{next} {gsub(/\r/,""); seq_len += length($0)} END{print seq_len}' "$fasta_file"
        fi
}

# Function to reverse complement a sequence
reverse_complement() {
        local sequence="$1"
        echo "$sequence" | rev | tr 'ACGTacgt' 'TGCAtgca'
}

# Function to optimize a single plasmid
optimize_plasmid() {
        local plasmid_file="$1"
        local plasmid_name
        plasmid_name=$(basename "$plasmid_file" .fasta)
        plasmid_name=$(basename "$plasmid_name" .fa)
        
        echo "Processing plasmid: $plasmid_name" >&2
        
        # Create temporary directory
        local temp_dir
        temp_dir="${OUTPUT_DIR}/temp_${plasmid_name}"
        mkdir -p "$temp_dir"
        
        # Create BLAST database for the plasmid
        makeblastdb -in "$plasmid_file" -dbtype nucl -out "${temp_dir}/plasmid_db" >/dev/null 2>&1
        
        # BLAST plasmid against reference
        blastn -query "$REFERENCE" -db "${temp_dir}/plasmid_db" \
                -out "${temp_dir}/blast_results.txt" \
                -outfmt "6 qseqid sseqid pident qcovs qstart qend sstart send slen qlen" \
                -evalue 1e-10 -num_threads "$THREADS"
        
        # Check if we have good alignments
        if [ ! -s "${temp_dir}/blast_results.txt" ]; then
                echo "  No significant alignments found" >&2
                # Try reverse complement
                echo "  Trying reverse complement..." >&2
                
                # Create reverse complement version
                local temp_rev="${temp_dir}/plasmid_rev.fasta"
                seqkit seq -r -p --seq-type DNA "$plasmid_file" > "$temp_rev" 2>/dev/null
                
                # BLAST reverse complement against reference
                makeblastdb -in "$temp_rev" -dbtype nucl -out "${temp_dir}/plasmid_rev_db" >/dev/null 2>&1
                blastn -query "$REFERENCE" -db "${temp_dir}/plasmid_rev_db" \
                        -out "${temp_dir}/blast_results_rev.txt" \
                        -outfmt "6 qseqid sseqid pident qcovs qstart qend sstart send slen qlen" \
                        -evalue 1e-10 -num_threads "$THREADS"
                
                if [ -s "${temp_dir}/blast_results_rev.txt" ]; then
                        echo "  Reverse complement shows better alignment" >&2
                        # Trim direct repeats before saving
                        local trimmed_rev="${temp_dir}/plasmid_rev_trimmed.fasta"
                        trim_direct_repeats "$temp_rev" "$trimmed_rev"
                        cp "$trimmed_rev" "${OUTPUT_DIR}/${plasmid_name}_optimized.fasta"
                        echo "  Optimized plasmid saved: ${plasmid_name}_optimized.fasta" >&2
                else
                        echo "  No good alignment even with reverse complement" >&2
                        cp "$plasmid_file" "${OUTPUT_DIR}/${plasmid_name}_unaligned.fasta"
                        echo "  Original plasmid saved: ${plasmid_name}_unaligned.fasta" >&2
                fi
        else
                # Analyze BLAST results
                local best_hit
                best_hit=$(sort -k3,3nr -k4,4nr "${temp_dir}/blast_results.txt" | head -n1)
                
                if [ -n "$best_hit" ]; then
                        local pident qcovs
                        pident=$(echo "$best_hit" | cut -f3)
                        qcovs=$(echo "$best_hit" | cut -f4)
                        
                        echo "  Best hit: ${pident}% identity, ${qcovs}% coverage" >&2
                        
                        # Process plasmid: check orientation, trim repeats, adjust start coordinates
                        local processed_file="${OUTPUT_DIR}/${plasmid_name}_processed.fasta"
                        local original_length=$(get_plasmid_length "$plasmid_file")
                        local repeat_bases_trimmed=0
                        
                        # Check orientation and create reverse complement if needed
                        echo "Checking orientation..." >&2
                        local temp_rev="${temp_dir}/plasmid_rev.fasta"
                        seqkit seq -r -p --seq-type DNA "$plasmid_file" > "$temp_rev" 2>/dev/null
                        
                        local better_orientation
                        better_orientation=$(get_better_orientation "$plasmid_file" "$temp_dir" "$temp_rev")
                        
                        local input_plasmid
                        if [ "$better_orientation" = "reverse" ]; then
                                if [ -f "$temp_rev" ]; then
                                        input_plasmid="$temp_rev"
                                        echo "  Using reverse complement" >&2
                                else
                                        echo "  WARNING: Reverse complement file not found, using original" >&2
                                        input_plasmid="$plasmid_file"
                                fi
                        else
                                input_plasmid="$plasmid_file"
                                echo "  Using original orientation" >&2
                        fi
                        
                        # Trim direct repeats
                        echo "Trimming direct repeats..." >&2
                        local trimmed_plasmid="${temp_dir}/plasmid_trimmed.fasta"
                        trim_direct_repeats "$input_plasmid" "$trimmed_plasmid"
                        
                        # Calculate repeat bases trimmed
                        local trimmed_length=$(get_plasmid_length "$trimmed_plasmid")
                        repeat_bases_trimmed=$((original_length - trimmed_length))
                        
                        # Find correct start coordinate
                        echo "Adjusting start coordinate..." >&2
                        find_correct_start "$trimmed_plasmid" "$processed_file"
                        
                        # Get final length
                        local final_length=$(get_plasmid_length "$processed_file")
                        
                        # Write statistics to file
                        echo -e "${plasmid_name}\t${original_length}\t${repeat_bases_trimmed}\t${final_length}\t${pident}\t${qcovs}" >> "$STATS_FILE"
                        
                        echo "Processed plasmid saved: ${plasmid_name}_processed.fasta" >&2
                fi
        fi
        
        # Clean up temporary files
        rm -rf "$temp_dir"
}

# Function to detect and trim direct repeat sequences at plasmid ends
trim_direct_repeats() {
        local plasmid_file="$1"
        local output_file="$2"
        
        echo "  Checking for direct repeat sequences at plasmid ends..." >&2
        
        # Read plasmid sequence
        local plasmid_seq
        plasmid_seq=$(awk 'BEGIN{seq=""} /^>/{next} {gsub(/\r/,""); seq = seq $0} END{print seq}' "$plasmid_file")
        local seq_len=${#plasmid_seq}
        
        if [ $seq_len -eq 0 ]; then
                echo "  Warning: Could not read plasmid sequence" >&2
                cp "$plasmid_file" "$output_file"
                return
        fi
        
        # Check for direct repeats at ends
        local max_repeat_len
        max_repeat_len=$((seq_len / 2))  # Can't be longer than half the sequence
        
        local best_repeat_len=0
        local best_repeat_start=0
        
        # Try different repeat lengths
        for ((repeat_len=10; repeat_len<=max_repeat_len; repeat_len++)); do
                local start_seq end_seq
                start_seq="${plasmid_seq:0:$repeat_len}"
                end_seq="${plasmid_seq: -$repeat_len}"
                
                if [ "$start_seq" = "$end_seq" ]; then
                        best_repeat_len=$repeat_len
                        best_repeat_start=$repeat_len
                        echo "  Found direct repeat of length $repeat_len at both ends" >&2
                fi
        done
        
        if [ $best_repeat_len -gt 0 ]; then
                echo "  Trimming $best_repeat_len bases from start (keeping end repeat)" >&2
            
                # Trim the repeat from the start, keep the end
                local trimmed_seq
                trimmed_seq="${plasmid_seq:$best_repeat_len}"
                
                # Write trimmed plasmid
                local plasmid_header
                plasmid_header=$(awk 'NR==1' "$plasmid_file")
                echo "$plasmid_header" > "$output_file"
                echo "$trimmed_seq" >> "$output_file"
                
                echo "  Direct repeat trimmed - new length: ${#trimmed_seq}" >&2
        else
                echo "  No significant direct repeats found" >&2
                cp "$plasmid_file" "$output_file"
        fi
}

# Function to test both orientations and find correct start coordinate
test_both_orientations() {
        local trimmed_plasmid="$1"
        local output_file="$2"
        
        # Find the correct start coordinate (orientation already determined)
        find_correct_start "$trimmed_plasmid" "$output_file"
}

# Function to find correct start coordinate by BLASTing reference start against plasmid
find_correct_start() {
        local plasmid_file="$1"
        local output_file="$2"
        
        echo "  Finding correct start coordinate..." >&2
        
        # Extract first 100 bp of reference
        local ref_start_100bp="${temp_dir}/ref_start_100bp.fasta"
        seqkit subseq -r 1:100 "$REFERENCE" > "$ref_start_100bp" 2>/dev/null
        
        # Create BLAST database for the plasmid
        makeblastdb -in "$plasmid_file" -dbtype nucl -out "${temp_dir}/plasmid_db" >/dev/null 2>&1
        
        # BLAST reference start against plasmid
        local blast_results="${temp_dir}/start_blast.txt"
        blastn -query "$ref_start_100bp" -db "${temp_dir}/plasmid_db" \
                -out "$blast_results" \
                -outfmt "6 qseqid sseqid pident qcovs qstart qend sstart send slen qlen" \
                -evalue 1e-10 -num_threads "$THREADS"
        
        if [ -s "$blast_results" ]; then
                # Get best hit
                local best_hit
                best_hit=$(sort -k3,3nr -k4,4nr "$blast_results" | head -n1)
                
                if [ -n "$best_hit" ]; then
                        # Parse coordinates: sstart is where reference start maps on plasmid
                        local plasmid_start
                        plasmid_start=$(echo "$best_hit" | cut -f7)
                        
                        echo "  Reference start (1-100) maps to plasmid position $plasmid_start" >&2
                        
                        # Adjust plasmid start coordinate
                        local plasmid_seq
                        plasmid_seq=$(awk 'BEGIN{seq=""} /^>/{next} {gsub(/\r/,""); seq = seq $0} END{print seq}' "$plasmid_file")
                        local seq_len=${#plasmid_seq}
                        
                        if [ $seq_len -gt 0 ] && [ "$plasmid_start" -gt 1 ]; then
                                # Calculate new start and end positions for circular permutation
                                local new_start=$((plasmid_start - 1))  # Convert to 0-based indexing
                                local new_end=$((seq_len - 1))
                                local new_seq
                                
                                # Create circular permutation: take sequence from new_start to end, then from beginning to new_start-1
                                new_seq="${plasmid_seq:$new_start}${plasmid_seq:0:$new_start}"
                                
                                # Save adjusted plasmid
                                local plasmid_header
                                plasmid_header=$(awk 'NR==1' "$plasmid_file")
                                echo "$plasmid_header" > "$output_file"
                                echo "$new_seq" >> "$output_file"
                                
                                echo "  Start coordinate adjusted to position $plasmid_start" >&2
                        else
                                echo "  No adjustment needed or invalid coordinates" >&2
                                cp "$plasmid_file" "$output_file"
                        fi
                else
                        echo "  No BLAST hits found for reference start" >&2
                        cp "$plasmid_file" "$output_file"
                fi
        else
                echo "  BLAST analysis failed" >&2
                cp "$plasmid_file" "$output_file"
        fi
        
        # Clean up temporary files
        rm -f "$ref_start_100bp" "${temp_dir}/plasmid_db"* "$blast_results"
}

# Function to get better orientation based on MUMmer coverage
get_better_orientation() {
        local plasmid_file="$1"
        local temp_dir="$2"
        local temp_rev="$3"  # Pass the reverse complement file path
        
        echo "  Checking for better orientation using MUMmer..." >&2
        
        # Run nucmer for both orientations
        nucmer --prefix="${temp_dir}/fwd" "$REFERENCE" "$plasmid_file" >/dev/null 2>&1
        nucmer --prefix="${temp_dir}/rev" "$REFERENCE" "$temp_rev" >/dev/null 2>&1
        
        # Get alignment coordinates and determine better orientation
        if [ -s "${temp_dir}/fwd.delta" ] && [ -s "${temp_dir}/rev.delta" ]; then
                local fwd_coords rev_coords
                fwd_coords=$(show-coords -r -c -l "${temp_dir}/fwd.delta" | tail -1)
                rev_coords=$(show-coords -r -c -l "${temp_dir}/rev.delta" | tail -1)
                
                if [ -n "$fwd_coords" ] && [ -n "$rev_coords" ]; then
                        # Parse coordinates: start2 < end2 = proper orientation
                        local fwd_start2 fwd_end2 rev_start2 rev_end2
                        fwd_start2=$(echo "$fwd_coords" | cut -d'|' -f2 | awk '{print $1}')
                        fwd_end2=$(echo "$fwd_coords" | cut -d'|' -f2 | awk '{print $2}')
                        rev_start2=$(echo "$rev_coords" | cut -d'|' -f2 | awk '{print $1}')
                        rev_end2=$(echo "$rev_coords" | cut -d'|' -f2 | awk '{print $2}')
                        
                        # Choose orientation with proper coordinate order (start2 < end2)
                        if [ "$fwd_start2" -lt "$fwd_end2" ] && [ "$rev_start2" -gt "$rev_end2" ]; then
                                echo "forward"
                        elif [ "$rev_start2" -lt "$rev_end2" ] && [ "$fwd_start2" -gt "$fwd_end2" ]; then
                                echo "reverse"
                        else
                                echo "forward" # Fallback to forward if no clear orientation
                        fi
                else
                        echo "forward" # Fallback to forward if no coordinates
                fi
        else
                echo "forward" # Fallback to forward if nucmer failed
        fi
        
        # Don't clean up temp files here - they might be needed by the caller
}

# Initialize statistics file
STATS_FILE="${OUTPUT_DIR}/processing_statistics.tsv"
echo -e "Sample_ID\tOriginal_length\tRepeat_bases_trimmed\tFinal_length\tBLAST_identity\tBLAST_coverage" > "$STATS_FILE"

# Find all plasmid files in the directory
echo "Scanning for plasmid files..." >&2
plasmid_files=()
shopt -s nullglob
for f in "$PLASMIDS_DIR"/*.fasta "$PLASMIDS_DIR"/*.fa "$PLASMIDS_DIR"/*.fas; do
        if [ -f "$f" ]; then
                plasmid_files+=("$f")
        fi
done
shopt -u nullglob

if [ ${#plasmid_files[@]} -eq 0 ]; then
        echo "No plasmid files found in: $PLASMIDS_DIR" >&2
        echo "Expected extensions: .fasta, .fa, .fas" >&2
        exit 1
fi

echo "Found ${#plasmid_files[@]} plasmid files" >&2

# Process each plasmid
total_plasmids=${#plasmid_files[@]}
current_plasmid=1

for plasmid_file in "${plasmid_files[@]}"; do
        echo "" >&2
        echo "==================================================================" >&2
        echo "Processing plasmid ${current_plasmid} of ${total_plasmids}: $(basename "$plasmid_file")" >&2
        echo "==================================================================" >&2
        
        if ! optimize_plasmid "$plasmid_file"; then
                echo "Failed to process plasmid: $(basename "$plasmid_file")" >&2
        fi
        
        current_plasmid=$((current_plasmid + 1))
done

echo "" >&2
echo "==================================================================" >&2
echo "OPTIMIZATION COMPLETE" >&2
echo "==================================================================" >&2
echo "Output directory: $OUTPUT_DIR" >&2
echo "Processed plasmids: ${total_plasmids}" >&2
echo "==================================================================" >&2
