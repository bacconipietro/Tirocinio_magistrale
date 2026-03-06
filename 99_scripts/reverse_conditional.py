import sys
import re

def reverse_complement(seq):
    """Reverse a DNA sequence"""
    return seq[::-1]

def parse_header(header):
    """Extract coordinates from header like >ID:c12628-11696"""
    # Look for pattern like c12628-11696 or similar coordinates
    match = re.search(r':c?(\d+)-(\d+)', header)
    if match:
        start = int(match.group(1))
        end = int(match.group(2))
        return start, end
    return None, None

def read_fasta(filename):
    """Read FASTA file and yield (header, sequence) tuples"""
    with open(filename, 'r') as f:
        header = None
        sequence = []
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    yield header, ''.join(sequence)
                header = line
                sequence = []
            else:
                sequence.append(line)
        
        if header:
            yield header, ''.join(sequence)

# Main script
if len(sys.argv) < 2:
    print("Usage: python3 reverse_conditional.py <input.fasta> [output.fasta]")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2] if len(sys.argv) > 2 else input_file.replace(".fasta", "_fixed.fasta")

with open(output_file, 'w') as out:
    for header, sequence in read_fasta(input_file):
        start, end = parse_header(header)
        
        # If start > end, reverse the sequence
        if start is not None and end is not None:
            if start > end:
                sequence = reverse_complement(sequence)
                print(f"Reversed: {header[:50]}... (start={start} > end={end})")
            else:
                print(f"Kept:     {header[:50]}... (start={start} < end={end})")
        
        # Write to output
        out.write(f"{header}\n")
        # Write sequence in 70-character lines for readability
        for i in range(0, len(sequence), 70):
            out.write(sequence[i:i+70] + '\n')

print(f"\n✓ Done! Output saved to: {output_file}")
