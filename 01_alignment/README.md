# Modify headers and align mitogenes with MAFFT  27-11-2025

### Modify headers, need to delete [..] because are not allowed in iqtree run.
```
for file in mt_genes_Mantodea/*.fasta; do
sed -i.old 's/_\[[^][]*\].*//g' "$file" 
done
```


## Align all fasta files with MAFFT

There are three main align options (localpair,genafpair,globalpair). In this situation i run --genafpair. The for cicle align all mtgenes files n one run, starting directory 01_mtgenes_Mantodea. **Before** starting the alignment process you must add **outgroups genes** for each file.
```
for file in mt_genes_Mantodea/*.fasta; do
base=(basename "$file" .fasta)
mafft --maxiterate 1000 --preservecase <align_option> "$file" >  "../01_alignments_Mantodea/${base}_alignment.fasta" 
done
```
align options: --localpair --genafpair --globalpair
