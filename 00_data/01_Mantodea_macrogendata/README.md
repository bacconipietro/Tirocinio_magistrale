# Trimmomatic e SPAdes 16-12-2025


-----

# Undestanding lcWGS outputs 30-01-2026

### Theory:

Your script is designed to run a "parameter sweep"—it tries to assemble the genome using increasing numbers of reads (50k, 100k, etc.) to find the minimum required for a good assembly.

It's data-driven: The assembler (MitoZ) looks for an overlap between the start and end of the sequence. If the coverage drops at the ends or the sequence is repetitive, it cannot confidently join them, so it outputs a linear sequence.

The Annotation consequence: Because the output file is technically linear, the Annotation step (which gives you the warning) treats it as a straight line. This is why your tRNAs are missing. If a gene sits exactly across the "break" point of the circle, the linear annotator cannot see it.



Check circularity in Ameass mt_assemble:

path sweep_single_ameass/runs/Ameass/20260123_114612_688230009/mt_sweep/reads_1000000

```bash
awk '/^>NODE_2($|[^0-9])/{print;p=1;next} /^>/{p=0} p' mtcandidate.fa > node_2.fasta
```
```bash
python3 -c "
import sys

# Change this to your actual filename if different
INPUT_FILE = 'node_2.fasta'

def check_circularity(file_path):
    seq = ''
    try:
        with open(file_path) as f:
            for line in f:
                line = line.strip()
                # Skip header lines, read only sequence
                if not line.startswith('>'):
                    seq += line
    except FileNotFoundError:
        print(f'Error: File {file_path} not found.')
        sys.exit(1)

    if not seq:
        print('Error: No sequence data found in file.')
        sys.exit(1)

    L = len(seq)
    print(f'Analyzing sequence (Length: {L} bp)')
    
    # Check for overlap at the ends (from 10bp up to 1000bp)
    overlap_len = 0
    # We loop up to 2000 just in case the overlap is large
    for i in range(10, 2000):
        # Check if the start (prefix) matches the end (suffix) of length i
        suffix = seq[-i:]
        if seq.startswith(suffix):
            overlap_len = i
    
    print(f'- Overlap detected: {overlap_len} bp')
    
    if overlap_len > 0:
        print(f'\nCONCLUSION: The genome IS CIRCULAR.') 
        print(f'To fix: You should trim {overlap_len} bp from the end of the sequence to avoid duplication.')
    else:
        print('\nCONCLUSION: No direct overlap found. It is likely linear or has a gap.')

if __name__ == '__main__':
    check_circularity(INPUT_FILE)
"
```


Output: CONCLUSION: No direct overlap found. It is likely linear or has a gap. 

```
Ameass: all outputs (mt.fasta,mt.gbf) 
1000000 OK      YES     YES     NODE_2  16189   106     313511  72.3694 100.00
Blastn: Select seq AY859585.1	Mantis religiosa internal transcribed spacer 1, partial sequence; 5.8S ribosomal RNA gene and internal transcribed spacer 2, complete sequence; and 28S ribosomal RNA gene, partial sequence	Mantis religiosa	6250	6522	(query cover=26%)	0.0	94.78%	4723	AY859585.1     NEGATIVE


Amedec: all outputs (mt.fasta,mt.gbf) 
2000000 OK      YES     YES     NODE_3  15664   297     808140  33.9074 99.73 
Blastn: Eremiaphila sp. mitochondrion, complete genome	Eremiaphila sp.	NA	2910531	9079	13645	(query cover=99%)	0.0	83.30%	15583	MG888444.1   POSITIVE


Aptapt: all outputs (mt.fasta,mt.gbf) 
2000000 OK      YES     YES     NODE_4  16427   217     612782  172.7433        100.00
Blastn: Select seq AY859585.1	Mantis religiosa internal transcribed spacer 1, partial sequence; 5.8S ribosomal RNA gene and internal transcribed spacer 2, complete sequence; and 28S ribosomal RNA gene, partial sequence	Mantis religiosa	6250	6522	(query cover=26%)	0.0	94.78%	4723	AY859585.1    		NEGATIVE


Amedum: runs failed, even 2 millions reads 
2000000 FAILED_VALIDATION_NO_HITS       NO      NO   No validated contigs detected in MitoZ summary.


Amepic: runs failed, 2000000 FAILED_MITOZ    NO   NO 


Amespa2: all outputs (mt.fasta,mt.gbf) 
1000000 OK      YES     YES     NODE_5  12293   166     359625  61.6058 99.47 
Blastn: Select seq OZ392354.1	Schoenobius gigantellus genome assembly, chromosome: 9	Schoenobius gigantellus	60.2	60.2	query cover=1%	0.010	81.33%	40936626	OZ392354.1


Pseben: all outputs (mt.fasta,mt.gbf) 
2000000 OK      YES     YES     NODE_1  19102   227     705455  36.1233 98.49
Blastn: Select seq OZ390435.1	Taeniapion urticarium genome assembly, chromosome: 8	Taeniapion urticarium	62.1	62.1	(query cover=0%)	0.004	100.00%	56091572	OZ390435.1


Rivbae: runs failed 2000000 FAILED_VALIDATION_NO_HITS   NO   NO  No validated contigs detected in MitoZ summary.
```
