#!/bin/bash
#
# Bash script to identify and flip sequences marked with complement() in headers
# Excludes tRNA sequences by default
# Uses mafft --adjustdirection for sequence flipping
#

set -e

# Configuration
VERBOSE=false
EXCLUDE_TRNA=true
FILE_PATTERN="*_fasta"
OUTPUT_DIR=""
INPUT_DIR=""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Print usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS] INPUT_DIR

Options:
  -o, --output DIR         Output directory (default: INPUT_DIR)
  -p, --pattern PATTERN    File pattern to match (default: *_fasta)
  -v, --verbose            Verbose output
  --no-exclude-trna        Do not skip tRNA sequences
  -h, --help               Show this help message

Examples:
  $0 /path/to/fasta/files
  $0 /path/to/fasta/files -o /path/to/output -v
  $0 /path/to/fasta/files -p "*.fasta" --no-exclude-trna

EOF
    exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -p|--pattern)
            FILE_PATTERN="$2"
            shift 2
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        --no-exclude-trna)
            EXCLUDE_TRNA=false
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            INPUT_DIR="$1"
            shift
            ;;
    esac
done

# Validate input
if [[ -z "$INPUT_DIR" ]]; then
    echo -e "${RED}Error: INPUT_DIR is required${NC}"
    usage
fi

if [[ ! -d "$INPUT_DIR" ]]; then
    echo -e "${RED}Error: $INPUT_DIR is not a directory${NC}"
    exit 1
fi

# Set output directory
if [[ -z "$OUTPUT_DIR" ]]; then
    OUTPUT_DIR="$INPUT_DIR"
fi

mkdir -p "$OUTPUT_DIR"

# Check for mafft
if ! command -v mafft &> /dev/null; then
    echo -e "${RED}Error: mafft is not installed${NC}"
    exit 1
fi

# Statistics
total_files=0
total_sequences=0
total_flipped=0
total_trna_skipped=0

# Function to check if header indicates complement
is_complement() {
    [[ $1 =~ \[location=complement\( ]]
}

# Function to check if header indicates tRNA
is_trna() {
    if [[ "$EXCLUDE_TRNA" == false ]]; then
        return 1
    fi
    # Check for common tRNA indicators
    if [[ $1 =~ _gt[0-9]+ ]] || [[ $1 =~ _tr[0-9]+ ]]; then
        return 0
    fi
    # Check for tRNA in gene field
    if [[ $1 =~ \[gene=tRNA ]]; then
        return 0
    fi
    return 1
}

# Function to extract sequence from FASTA file
extract_sequence_only() {
    grep -v "^>" "$1"
}

# Function to extract header from FASTA file
extract_header_only() {
    grep "^>" "$1" | head -1
}

# Function to process a single file
process_file() {
    local input_file=$1
    local output_file=$2
    local filename=$(basename "$input_file")
    
    if [[ "$VERBOSE" == true ]]; then
        echo -e "\n${GREEN}Processing: $filename${NC}"
    fi
    
    local file_total=0
    local file_flipped=0
    local file_trna_skipped=0
    local temp_output=$(mktemp)
    
    local current_header=""
    local current_seq=""
    local in_seq=false
    
    # Process each sequence in the file
    while IFS= read -r line; do
        if [[ $line == ">"* ]]; then
            # New sequence - process previous one if exists
            if [[ -n "$current_header" ]]; then
                file_total=$((file_total + 1))
                
                local has_complement=false
                local is_trna=false
                
                is_complement "$current_header" && has_complement=true
                is_trna "$current_header" && is_trna=true
                
                if [[ "$has_complement" == true && "$is_trna" == true ]]; then
                    if [[ "$VERBOSE" == true ]]; then
                        echo -e "  ${YELLOW}SKIP (tRNA + complement)${NC}: ${current_header:0:60}..."
                    fi
                    file_trna_skipped=$((file_trna_skipped + 1))
                    echo "$current_header" >> "$temp_output"
                    echo "$current_seq" >> "$temp_output"
                elif [[ "$has_complement" == true ]]; then
                    if [[ "$VERBOSE" == true ]]; then
                        echo -e "  ${GREEN}FLIP${NC}: ${current_header:0:60}..."
                    fi
                    
                    # Create temporary FASTA for mafft
                    local temp_seq=$(mktemp)
                    echo "$current_header" > "$temp_seq"
                    echo "$current_seq" >> "$temp_seq"
                    
                    # Run mafft and capture flipped sequence
                    if flipped_seq=$(mafft --adjustdirection "$temp_seq" 2>/dev/null | grep -v "^>"); then
                        echo "$current_header" >> "$temp_output"
                        echo "$flipped_seq" >> "$temp_output"
                        file_flipped=$((file_flipped + 1))
                    else
                        echo -e "${YELLOW}WARNING: Failed to flip, keeping original${NC}: ${current_header:0:60}..." >&2
                        echo "$current_header" >> "$temp_output"
                        echo "$current_seq" >> "$temp_output"
                    fi
                    rm -f "$temp_seq"
                else
                    echo "$current_header" >> "$temp_output"
                    echo "$current_seq" >> "$temp_output"
                fi
            fi
            
            # Start new sequence
            current_header="$line"
            current_seq=""
            in_seq=true
        else
            # Continuation of sequence
            current_seq="${current_seq}${line}"
        fi
    done < "$input_file"
    
    # Process last sequence
    if [[ -n "$current_header" ]]; then
        file_total=$((file_total + 1))
        
        local has_complement=false
        local is_trna=false
        
        is_complement "$current_header" && has_complement=true
        is_trna "$current_header" && is_trna=true
        
        if [[ "$has_complement" == true && "$is_trna" == true ]]; then
            if [[ "$VERBOSE" == true ]]; then
                echo -e "  ${YELLOW}SKIP (tRNA + complement)${NC}: ${current_header:0:60}..."
            fi
            file_trna_skipped=$((file_trna_skipped + 1))
            echo "$current_header" >> "$temp_output"
            echo "$current_seq" >> "$temp_output"
        elif [[ "$has_complement" == true ]]; then
            if [[ "$VERBOSE" == true ]]; then
                echo -e "  ${GREEN}FLIP${NC}: ${current_header:0:60}..."
            fi
            
            # Create temporary FASTA for mafft
            local temp_seq=$(mktemp)
            echo "$current_header" > "$temp_seq"
            echo "$current_seq" >> "$temp_seq"
            
            # Run mafft and capture flipped sequence
            if flipped_seq=$(mafft --adjustdirection "$temp_seq" 2>/dev/null | grep -v "^>"); then
                echo "$current_header" >> "$temp_output"
                echo "$flipped_seq" >> "$temp_output"
                file_flipped=$((file_flipped + 1))
            else
                echo -e "${YELLOW}WARNING: Failed to flip, keeping original${NC}: ${current_header:0:60}..." >&2
                echo "$current_header" >> "$temp_output"
                echo "$current_seq" >> "$temp_output"
            fi
            rm -f "$temp_seq"
        else
            echo "$current_header" >> "$temp_output"
            echo "$current_seq" >> "$temp_output"
        fi
    fi
    
    # Move temp output to final output
    mv "$temp_output" "$output_file"
    
    if [[ "$VERBOSE" == true ]]; then
        echo -e "  Summary: Total=$file_total, Flipped=$file_flipped, Skipped tRNA=$file_trna_skipped"
        echo -e "  Output: $output_file"
    fi
    
    # Update totals
    total_sequences=$((total_sequences + file_total))
    total_flipped=$((total_flipped + file_flipped))
    total_trna_skipped=$((total_trna_skipped + file_trna_skipped))
}

# Find and process all files
cd "$INPUT_DIR"
for file in $FILE_PATTERN; do
    if [[ -f "$file" ]]; then
        total_files=$((total_files + 1))
        output_file="$OUTPUT_DIR/$file"
        process_file "$file" "$output_file"
    fi
done

# Print summary
echo -e "\n${GREEN}$(printf '%0.s=' {1..70})${NC}"
echo -e "${GREEN}SUMMARY${NC}"
echo -e "${GREEN}$(printf '%0.s=' {1..70})${NC}"
echo "Files processed:          $total_files"
echo "Total sequences:          $total_sequences"
echo "Sequences flipped:        $total_flipped"
echo "Sequences skipped (tRNA): $total_trna_skipped"
echo -e "Output directory:         $OUTPUT_DIR"
echo -e "${GREEN}$(printf '%0.s=' {1..70})${NC}"
