#!/bin/bash

################################################################################
# Script: fix_reverse_complements.sh
# Purpose: Find FASTA files containing reverse complement sequences (marked by :c)
#          and flip them using MAFFT's --adjustdirection option
# Usage: ./fix_reverse_complements.sh [search_directory OR fasta_files...]
#        ./fix_reverse_complements.sh old/
#        ./fix_reverse_complements.sh old/*.fasta
#        ./fix_reverse_complements.sh file1.fasta file2.fasta file3.fasta
################################################################################

set -euo pipefail

# Configuration
OUTPUT_SUFFIX="_rc_fixed"
LOG_FILE="reverse_complement_fixes.log"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Initialize log
echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting reverse complement detection and fixing" | tee "$LOG_FILE"

# Check if MAFFT is installed
if ! command -v mafft &> /dev/null; then
    echo -e "${RED}ERROR: mafft is not installed or not in PATH${NC}" | tee -a "$LOG_FILE"
    exit 1
fi

echo -e "${BLUE}MAFFT found: $(which mafft)${NC}" | tee -a "$LOG_FILE"
echo "---" | tee -a "$LOG_FILE"

# Counter variables
total_files=0
files_with_rc=0
files_processed=0
files_failed=0

# Function to process a single FASTA file
process_fasta_file() {
    local fasta_file="$1"
    
    if [[ ! -f "$fasta_file" ]]; then
        echo -e "${RED}  ✗ File not found: $fasta_file${NC}" | tee -a "$LOG_FILE"
        files_failed=$((files_failed + 1))
        return 1
    fi
    
    total_files=$((total_files + 1))
    
    # Check if file contains reverse complement markers (:c)
    if grep -q ":c" "$fasta_file"; then
        files_with_rc=$((files_with_rc + 1))
        
        echo -e "${YELLOW}[$(date '+%H:%M:%S')] Processing: $fasta_file${NC}" | tee -a "$LOG_FILE"
        
        # Create output filename
        base_name="${fasta_file%.*}"
        extension="${fasta_file##*.}"
        output_file="${base_name}${OUTPUT_SUFFIX}.${extension}"
        
        # Check if output file already exists
        if [[ -f "$output_file" ]]; then
            echo -e "${YELLOW}  ⚠ Output file already exists: $output_file${NC}" | tee -a "$LOG_FILE"
            echo "  Backing up to ${output_file}.bak" | tee -a "$LOG_FILE"
            mv "$output_file" "${output_file}.bak"
        fi
        
        # Run MAFFT with --adjustdirection
        if mafft --adjustdirection "$fasta_file" > "$output_file" 2>> "$LOG_FILE"; then
            files_processed=$((files_processed + 1))
            
            # Count sequences in output
            seq_count=$(grep "^>" "$output_file" | wc -l)
            echo -e "${GREEN}  ✓ Success: Created $output_file (${seq_count} sequences)${NC}" | tee -a "$LOG_FILE"
            
        else
            files_failed=$((files_failed + 1))
            echo -e "${RED}  ✗ Failed: Could not process $fasta_file${NC}" | tee -a "$LOG_FILE"
            # Remove incomplete output file
            [[ -f "$output_file" ]] && rm "$output_file"
            return 1
        fi
    else
        echo -e "${BLUE}[$(date '+%H:%M:%S')] Skipping (no :c markers): $fasta_file${NC}" | tee -a "$LOG_FILE"
    fi
}

# Check arguments
if [[ $# -eq 0 ]]; then
    # No arguments, search current directory
    echo "Search directory: ." | tee -a "$LOG_FILE"
    
    while IFS= read -r fasta_file; do
        if [[ -n "$fasta_file" ]]; then
            process_fasta_file "$fasta_file"
        fi
    done < <(find "." -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.faa" -o -name "*.fna" \) 2>/dev/null)
    
elif [[ -d "$1" ]]; then
    # First argument is a directory
    SEARCH_DIR="$1"
    echo "Search directory: $SEARCH_DIR" | tee -a "$LOG_FILE"
    
    while IFS= read -r fasta_file; do
        if [[ -n "$fasta_file" ]]; then
            process_fasta_file "$fasta_file"
        fi
    done < <(find "$SEARCH_DIR" -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.faa" -o -name "*.fna" \) 2>/dev/null)
    
else
    # Arguments are files
    echo "Processing specified files:" | tee -a "$LOG_FILE"
    for fasta_file in "$@"; do
        process_fasta_file "$fasta_file"
    done
fi

# Summary
echo "" | tee -a "$LOG_FILE"
echo "---" | tee -a "$LOG_FILE"
echo -e "${BLUE}SUMMARY:${NC}" | tee -a "$LOG_FILE"
echo "Total FASTA files processed: $total_files" | tee -a "$LOG_FILE"
echo -e "${YELLOW}Files with reverse complement markers (:c): $files_with_rc${NC}" | tee -a "$LOG_FILE"
echo -e "${GREEN}Files successfully fixed: $files_processed${NC}" | tee -a "$LOG_FILE"
echo -e "${RED}Files failed: $files_failed${NC}" | tee -a "$LOG_FILE"
echo "Output files saved with suffix: ${OUTPUT_SUFFIX}" | tee -a "$LOG_FILE"
echo "Log file: $LOG_FILE" | tee -a "$LOG_FILE"
echo "$(date '+%Y-%m-%d %H:%M:%S') - Completed" | tee -a "$LOG_FILE"

# Exit with success if all processing succeeded
[[ $files_failed -eq 0 ]] && exit 0 || exit 1
