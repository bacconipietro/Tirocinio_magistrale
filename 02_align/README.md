# Alignment and Trimming pipeline
## Aligning with MAFFT v7.525
Align nucleotide sequence 10-03-2026

```bash
mkdir 00_genafpair
mkdir 01_globalpair
mkdir 02_localpair
for file in /DATABIG/pietrobacconi/ncbi_datasets/refseq_mitogenomes/01_mtGENES/mt_genes/*.fasta; do
 base=$(basename "$file" .fasta)
 mafft --adjustdirection --maxiterate 1000 --preservecase --genafpair "$file" >  "00_genafpair/${base}_aligned.fasta" 
 mafft --adjustdirection --maxiterate 1000 --preservecase --globalpair "$file" >  "01_globalpair/${base}_aligned.fasta" 
 mafft --adjustdirection --maxiterate 1000 --preservecase --localpair "$file" >  "02_localpair/${base}_aligned.fasta" 
done
```
Align AA sequence

```bash
mkdir 00_genafpair
mkdir 01_globalpair
mkdir 02_localpair
for file in /DATABIG/pietrobacconi/ncbi_datasets/refseq_mitogenomes/01_mtGENES/AA/*.fasta; do
 base=$(basename "$file" .fasta)
 mafft --adjustdirection --maxiterate 1000 --preservecase --genafpair "$file" >  "00_genafpair/${base}_aligned.fasta" 
 mafft --adjustdirection --maxiterate 1000 --preservecase --globalpair "$file" >  "01_globalpair/${base}_aligned.fasta" 
 mafft --adjustdirection --maxiterate 1000 --preservecase --localpair "$file" >  "02_localpair/${base}_aligned.fasta" 
done
```
-----

## Trimming with transeq (EMBOSS:6.6.0.0)

```bash
for file in ../../../02_Alignments/00_nucleotide/<align_mode>/*_aligned.fasta; do
 trimal -in "$file" -out ${file/_aligned.fasta/_trimmed.fasta} -automated1 -keepheader -fasta 
done
```

```bash
for file in ../../../02_Alignments/01_AA/<align_mode>/*_aligned.fasta; do
 trimal -in "$file" -out ${file/_aligned.fasta/_trimmed.fasta} -automated1 -keepheader -fasta 
done
```

After trimming we need to remove `_R_` pattern which is pasted by MAFFT flag when it reverses a sequence.
```bash
for ile in *_trimmed.fasta; do
sed -i 's/_R_//' "$file"
done
```
