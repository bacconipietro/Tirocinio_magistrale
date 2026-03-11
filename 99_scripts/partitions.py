import re
import sys

input_file = sys.argv[1]
output_file = sys.argv[1].rsplit(".", 1)[0] + "_codon." + sys.argv[1].rsplit(".", 1)[1]

with open(input_file) as f:
    lines = f.readlines()

out = []
counter = 1

for line in lines:
    m = re.match(r'\s*charset\s+p\d+_(\w+?)_trimmed\s*=\s*(\d+)-(\d+);', line)
    if m:
        gene = m.group(1)
        start = int(m.group(2))
        end = int(m.group(3))

        if gene.lower().startswith("rrn"):
            # Ribosomal gene: single line, no codon partitioning
            out.append(f"\tcharset p{counter}_{gene} = {start}-{end};")
            counter += 1
        else:
            out.append(f"\tcharset p{counter}_{gene} = {start}-{end}\\3;")
            out.append(f"\tcharset p{counter+1}_{gene} = {start+1}-{end}\\3;")
            out.append(f"\tcharset p{counter+2}_{gene} = {start+2}-{end}\\3;")
            counter += 3
    else:
        out.append(line.rstrip())

with open(output_file, "w") as f:
    f.write("\n".join(out) + "\n")

print(f"Done! Output written to {output_file}")


