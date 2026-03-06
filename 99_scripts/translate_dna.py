from Bio import SeqIO
import sys

# Check if input file is provided
if len(sys.argv) < 2:
    print("Usage: python3 translate_dna.py <input.fasta> [output.fasta]")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2] if len(sys.argv) > 2 else input_file.replace(".fasta", "_AA.fasta")

# Read and translate
with open(output_file, 'w') as out:
    for record in SeqIO.parse(input_file, "fasta"):
        protein_seq = str(record.seq.translate())
        out.write(f">{record.id}\n{protein_seq}\n")

print(f"✓ Translated: {input_file} → {output_file}")
