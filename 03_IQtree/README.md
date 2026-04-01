# IQtree pipeline

-----

# AMAS 11-03-2026
Make a concat file with trimmed alignements. First of all we need to remove _R_ pattern from headers, it's very important to do at first because concat command nedds to have same headers. 

```bash
for f in 03_Trimming/00_nucleotide/00_genafpair/*.fasta; do
    sed 's/_R_//' "$f" > "04_IQtree/00_concat/00_genafpair/$(basename "$f")"
done
```


Then run AMAS to output partitions.nex, giving in input only the alignment (genafpair)

```bash
AMAS.py concat -p nu_genaf_partitions.nex -t nu_genaf_concat.fasta -u fasta -y nexus -i *.fasta -f fasta -d dna 
```

After all we run a homemade script to define partitions for each gene, the following model selection will find a model for each partition. 

```bash
python3 /home/.../partitions.py nu_gneaf_partitions.txt 
```
