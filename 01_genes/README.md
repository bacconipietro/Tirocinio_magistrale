# Mito Genes giant grep 10-03-2026

### Grep
New folders structure with all simplify directories linked in 

+ downloads
+ assembly
+ outgroups 

Total number of species attended = 189

```bash
mkdir -p mt_genes
#the flag -L is specific for links)
for file in $(find -L . -name "*.fasta"); do     
    awk '
        /^>/ {
            p = 0
            low = tolower($0)

            if (low ~ /gene=(cox1|coxi([^a-z]|$)|coi([^a-z]|$))/)       { out="mt_genes/COX1.fasta"; p=1 }
            else if (low ~ /gene=(cox2|coxii([^a-z]|$)|coii([^a-z]|$))/)       { out="mt_genes/COX2.fasta"; p=1 }
            else if (low ~ /gene=(cox3|coxiii|coiii)/)      { out="mt_genes/COX3.fasta"; p=1 }
            else if (low ~ /gene=(cob|cytb|cyt_b)/)        { out="mt_genes/CYTB.fasta"; p=1 }
            else if (low ~ /gene=(atp6)/)            { out="mt_genes/ATP6.fasta"; p=1 }
            else if (low ~ /gene=(atp8)/)            { out="mt_genes/ATP8.fasta"; p=1 }
            else if (low ~ /gene=(nd1|nad1)/)        { out="mt_genes/ND1.fasta"; p=1 }
            else if (low ~ /gene=(nd2|nad2)/)        { out="mt_genes/ND2.fasta"; p=1 }
            else if (low ~ /gene=(nd3|nad3)/)        { out="mt_genes/ND3.fasta"; p=1 }
            else if (low ~ /gene=(nd4([^a-z]|$)|nad4([^a-z]|$))/)             { out="mt_genes/ND4.fasta"; p=1 }
            else if (low ~ /gene=(nd4l|nad4l)/) { out="mt_genes/ND4L.fasta"; p=1 }
            else if (low ~ /gene=(nd5|nad5)/)        { out="mt_genes/ND5.fasta"; p=1 }
            else if (low ~ /gene=(nd6|nad6)/)        { out="mt_genes/ND6.fasta"; p=1 }
            else if (low ~ /(gene|locus_tag)=(rrn12|rrns|rrnS|s-rrna|.*_mgr01|.*_gr01)/) { out="mt_genes/rrns.fasta"; p=1 }
            else if (low ~ /(gene|locus_tag)=(rrn16|rrnl|rrnL|l-rrna|.*_mgr02|.*_gr02)/) { out="mt_genes/rrnl.fasta"; p=1 }
            else { out="mt_genes/UNMATCHED_genes.fasta"; p=1 }
        }

        p { print >> out }
    ' "$file"
done
```

### Translate in AA

```bash
conda install bioconda::emboss

transeq -sequence ATP6.fasta -outseq ATP6_AA.fasta -frame 1 -table 5 -trim -clean 

for file in mt_genes/*.fasta; do
    base=$(basename "$file" .fasta)
    transeq -sequence "$file" -outseq "AA/${base}_AA.fasta" -frame 1 -table 5 -trim -clean
done
```

### Details to go on with pipeline

```bash
# Check genes count
for file in mt_genes/*.fasta; do

    count=$(grep "^>" "$file" | wc -l)
    echo "Tips number for $file is $count"
done 

# Check possible duplicates
grep "^>" your_file.fasta | sort | uniq -d


#delete _[]_ because are not allowed in IQtree
for file in *.fasta; do
sed -i.old 's/_\[[^][]*\].*//g' "$file" 
done
```
