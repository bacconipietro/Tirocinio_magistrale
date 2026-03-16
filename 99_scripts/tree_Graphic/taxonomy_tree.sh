#!/bin/bash

# Script to annotate a Newick tree file with family information from a TSV taxonomy file
# Usage: ./annotate_tree.sh input.nwk taxonomy.tsv output.nwk

if [ $# -ne 3 ]; then
    echo "Usage: $0 <input.nwk> <taxonomy.tsv> <output.nwk>"
    exit 1
fi

INPUT_NWK="$1"
TAX_TSV="$2"
OUTPUT_NWK="$3"

# Create a temporary awk script that builds a lookup table and performs substitutions
awk -v nwk_file="$INPUT_NWK" '
BEGIN {
    # Read the Newick file
    while ((getline line < nwk_file) > 0) {
        newick = line
    }
    close(nwk_file)
}

NR == 1 {
    # Skip header line
    next
}

{
    # Build associative array: accession_code -> family|subfamily
    # Note: order is ACCESSION_CODE (as it appears in contree file)
    accession_code = $2 "_" $1
    family = $3
    subfamily = $4
    tax[accession_code] = family "|" subfamily
}

END {
    # Now process the newick string and replace each code_accession with code_accession_family
    result = newick
    
    # Iterate through all entries in the tax array
    for (key in tax) {
        # Use gsub to replace all occurrences of the pattern
        # Pattern: the key followed by either ":" or "," or ")" or ";" (end of branch label)
        pattern = key "([,):]|;)"
        replacement = key "_" tax[key] "\\1"
        
        # We need to use a more robust approach with character-by-character processing
        # or use a sed-like approach
    }
    
    # Better approach: use character-by-character parsing
    output = ""
    i = 1
    while (i <= length(result)) {
        # Try to match an accession_code at current position
        matched = 0
        for (key in tax) {
            if (substr(result, i, length(key)) == key) {
                # Check that it'\''s a complete match (followed by delimiter)
                next_char = substr(result, i + length(key), 1)
                if (next_char == ":" || next_char == "," || next_char == ")" || next_char == ";") {
                    # Split family|subfamily and add both
                    split(tax[key], parts, "|")
                    output = output key "_" parts[1] "_" parts[2]
                    i += length(key)
                    matched = 1
                    break
                }
            }
        }
        if (!matched) {
            output = output substr(result, i, 1)
            i++
        }
    }
    
    print output
}
' "$TAX_TSV" > "$OUTPUT_NWK"

echo "Tree annotation complete!"
echo "Input: $INPUT_NWK"
echo "Taxonomy: $TAX_TSV"
echo "Output: $OUTPUT_NWK"
