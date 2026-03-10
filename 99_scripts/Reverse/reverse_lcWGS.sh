#!/bin/bash

################################################################################
# Script: flip_colon_minus_strand_sequences.sh
# Purpose: Find FASTA headers containing ":-" pattern (indicating reverse strand)
#          and flip the corresponding sequences to forward orientation
# Usage: ./flip_colon_minus_strand_sequences.sh [fasta_files...]
#        ./flip_colon_minus_strand_sequences.sh *.fasta
#        ./flip_colon_minus_strand_sequences.sh file1.fasta file2.fasta
################################################################################

set -euo pipefail

# Configuration
OUTPUT_SUFFIX="_flipped"
LOG_FILE="sequence_flips.log"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Initialize log
echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting colon-minus strand sequence flipping" | tee "$LOG_FILE"
echo "---" | tee -a "$LOG_FILE"

# Counter variables
total_files=0
files_with_minus=0
sequences_flipped=0
files_failed=0

# Function to reverse complement a DNA sequence
reverse_complement() {
    local seq="$1"
    # Reverse the sequence and complement each base
    echo "$seq" | rev | tr 'ACGTacgt' 'TGCAtgca'
}

# Function to process a single FASTA file
process_fasta_file() {
    local input_file="$1"
    local output_file="${input_file%.*}${OUTPUT_SUFFIX}.${input_file##*.}"
    local temp_file="${output_file}.tmp"
    
    if [[ ! -f "$input_file" ]]; then
        echo -e "${RED}  ✗ File not found: $input_file${NC}" | tee -a "$LOG_FILE"
        files_failed=$((files_failed + 1))
        return 1
    fi
    
    total_files=$((total_files + 1))
    
    # Check if file contains colon-minus strand markers (:-) in headers
    if ! grep -q ":-" "$input_file"; then
        echo -e "${BLUE}[$(date '+%H:%M:%S')] Skipping (no :- markers): $input_file${NC}" | tee -a "$LOG_FILE"
        return 0
    fi
    
    files_with_minus=$((files_with_minus + 1))
    echo -e "${YELLOW}[$(date '+%H:%M:%S')] Processing: $input_file${NC}" | tee -a "$LOG_FILE"
    
    # Check if output file already exists
    if [[ -f "$output_file" ]]; then
        echo -e "${YELLOW}  ⚠ Output file already exists: $output_file${NC}" | tee -a "$LOG_FILE"
        echo "  Backing up to ${output_file}.bak" | tee -a "$LOG_FILE"
        mv "$output_file" "${output_file}.bak"
    fi
    
    # Process the FASTA file
    local current_header=""
    local current_sequence=""
    local local_flipped=0
    
    {
        while IFS= read -r line; do
            if [[ "$line" =~ ^">" ]]; then
                # If we have a previous sequence, write it out
                if [[ -n "$current_header" ]]; then
                    if [[ "$current_header" =~ \:- ]]; then
                        # This header has colon-minus strand marker - flip the sequence
                        flipped_seq=$(reverse_complement "$current_sequence")
                        echo "$current_header" >> "$temp_file"
                        echo "$flipped_seq" >> "$temp_file"
                        local_flipped=$((local_flipped + 1))
                    else
                        # Forward strand - write as is
                        echo "$current_header" >> "$temp_file"
                        echo "$current_sequence" >> "$temp_file"
                    fi
                fi
                
                # Start new sequence
                current_header="$line"
                current_sequence=""
            else
                # Append to current sequence (remove whitespace)
                current_sequence="${current_sequence}${line//[[:space:]]/}"
            fi
        done
        
        # Don't forget the last sequence
        if [[ -n "$current_header" ]]; then
            if [[ "$current_header" =~ \:- ]]; then
                flipped_seq=$(reverse_complement "$current_sequence")
                echo "$current_header" >> "$temp_file"
                echo "$flipped_seq" >> "$temp_file"
                local_flipped=$((local_flipped + 1))
            else
                echo "$current_header" >> "$temp_file"
                echo "$current_sequence" >> "$temp_file"
            fi
        fi
    } < "$input_file"
    
    # Verify temp file was created and has content
    if [[ ! -f "$temp_file" ]] || [[ ! -s "$temp_file" ]]; then
        echo -e "${RED}  ✗ Failed: Could not create output file${NC}" | tee -a "$LOG_FILE"
        files_failed=$((files_failed + 1))
        [[ -f "$temp_file" ]] && rm "$temp_file"
        return 1
    fi
    
    # Move temp file to final output
    mv "$temp_file" "$output_file"
    sequences_flipped=$((sequences_flipped + local_flipped))
    
    # Count total sequences
    seq_count=$(grep "^>" "$output_file" | wc -l)
    echo -e "${GREEN}  ✓ Success: Created $output_file${NC}" | tee -a "$LOG_FILE"
    echo "    - Sequences flipped: $local_flipped" | tee -a "$LOG_FILE"
    echo "    - Total sequences: $seq_count" | tee -a "$LOG_FILE"
}

# Check if files are provided
if [[ $# -eq 0 ]]; then
    echo -e "${RED}ERROR: No files specified${NC}" | tee -a "$LOG_FILE"
    echo "Usage: $0 file1.fasta file2.fasta ..." | tee -a "$LOG_FILE"
    exit 1
fi

# Process each file
echo "Processing specified files:" | tee -a "$LOG_FILE"
for fasta_file in "$@"; do
    process_fasta_file "$fasta_file"
done

# Summary
echo "" | tee -a "$LOG_FILE"
echo "---" | tee -a "$LOG_FILE"
echo -e "${BLUE}SUMMARY:${NC}" | tee -a "$LOG_FILE"
echo "Total FASTA files processed: $total_files" | tee -a "$LOG_FILE"
echo -e "${YELLOW}Files with colon-minus strand markers (:-): $files_with_minus${NC}" | tee -a "$LOG_FILE"
echo -e "${GREEN}Total sequences flipped: $sequences_flipped${NC}" | tee -a "$LOG_FILE"
echo -e "${RED}Files failed: $files_failed${NC}" | tee -a "$LOG_FILE"
echo "Output files saved with suffix: ${OUTPUT_SUFFIX}" | tee -a "$LOG_FILE"
echo "Log file: $LOG_FILE" | tee -a "$LOG_FILE"
echo "$(date '+%Y-%m-%d %H:%M:%S') - Completed" | tee -a "$LOG_FILE"

# Exit with success if all processing succeeded
[[ $files_failed -eq 0 ]] && exit 0 || exit 1
