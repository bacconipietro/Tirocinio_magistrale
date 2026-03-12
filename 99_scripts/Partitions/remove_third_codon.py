import re
import sys
import os

def remove_third_codon_partitions(input_file):
    with open(input_file, 'r') as f:
        content = f.read()

    lines = content.split('\n')
    output_lines = []
    removed_count = 0
    kept_count = 0

    # First pass: collect gene starts for all codon-partitioned genes
    gene_starts = {}
    for line in lines:
        charset_match = re.match(r'\s*charset\s+p\d+(_\w+)\s*=\s*(\d+)-\d+\\3;', line, re.IGNORECASE)
        if charset_match:
            gene_tag = charset_match.group(1)
            start = int(charset_match.group(2))
            if gene_tag not in gene_starts:
                gene_starts[gene_tag] = start
            else:
                gene_starts[gene_tag] = min(gene_starts[gene_tag], start)

    # Second pass: remove 3rd codon position lines, keep everything else untouched
    for line in lines:
        charset_match = re.match(r'(\s*charset\s+)(p\d+)(_\w+)\s*=\s*(\S+);', line, re.IGNORECASE)

        if charset_match:
            gene_tag = charset_match.group(3)
            position = charset_match.group(4)

            is_codon = '\\3' in position

            if is_codon:
                range_match = re.match(r'(\d+)-\d+\\3', position)
                if range_match:
                    start = int(range_match.group(1))
                    gene_start = gene_starts.get(gene_tag, start)
                    offset = (start - gene_start) % 3  # 0=1st, 1=2nd, 2=3rd

                    if offset == 2:  # 3rd codon position → remove
                        removed_count += 1
                        continue  # Skip this line, keep everything else untouched

            # Keep line completely untouched
            output_lines.append(line)
            kept_count += 1
        else:
            output_lines.append(line)

    # Build output filename
    base, ext = os.path.splitext(input_file)
    output_file = base + '_W3' + ext

    with open(output_file, 'w') as f:
        f.write('\n'.join(output_lines))

    print(f"Done! Kept {kept_count} partitions, removed {removed_count} third-codon partitions.")
    print(f"Output written to: {output_file}")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python remove_third_codon.py input_file.nex")
        sys.exit(1)

    input_file = sys.argv[1]
    if not os.path.exists(input_file):
        print(f"Error: File '{input_file}' not found.")
        sys.exit(1)

    remove_third_codon_partitions(input_file)
