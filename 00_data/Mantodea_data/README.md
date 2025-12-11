# Download NCBI dataset 18-11-2025

```
mkdir -p downloads downloads_PCGs downloads_failed downloads_txt
set -uo pipefail
 
cut -f2,3 Mantodea_NCBIdataset.tsv | tail -n +2 - | while IFS=$'\t' read species acc;
do
    [ -z "$acc" ] && continue
    echo "Processing ${acc} for ${species}..."
    if ! esearch -db nucleotide -query "${acc}" </dev/null | efetch -format gene_fasta > "${acc}_${species}.gene_fasta"; then
       echo "${acc}_${species}" >> downloads_txt/failed_downloads.txt
       echo "Download "${acc}_${species}" failed"
       continue
    fi
    [ -s "${acc}_${species}.gene_fasta" ] || {
    echo "${acc}_${species}.gene_fasta is empty"
    echo "${acc}_${species}" >> downloads_txt/failed_downloads.txt
    rm "${acc}_${species}.gene_fasta"
    continue
    } 
       count=$( (grep -c "^>" "${acc}_${species}.gene_fasta" 2>/dev/null) || echo 0 )    
    if  [ "$count" -le 13 ]; then
         sed -i '/^>/ s/ /_/g' ${acc}_${species}.gene_fasta
         echo "${acc}_${species}" >> downloads_txt/downloads_only_PCGs.txt
         echo "${acc}_${species} has only $count sequences"
         mv "${acc}_${species}.gene_fasta" downloads_PCGs/
       else
       sed -i '/^>/ s/ /_/g' ${acc}_${species}.gene_fasta
       echo "${acc}_${species} downloaded correctly with $count sequences"
       echo "${acc}_${species}" >> downloads_txt/correct_downloads.txt
       mv "${acc}_${species}.gene_fasta" downloads/ 
    fi
done
```
# Simplify headers 19-11-2025

### Upload di Mantodea_rrnl.fasta Mantodea_rrns.fasta
```
mkdir Download_MantodeaNCBIdataset/downloads_rRNAs
scp Desktop/data/rRNAs_downloads/*.fasta  STUDENTI^pietro.bacconi@137.204.142.152:/DATABIG/pietrobacconi/ncbi_datasets/refseq_mitogenomes/Download_MantodeaNCBIdataset/downloads_rRNAs/
```

### Upload di AA001_Ameand.fasta AA003_Ameser.fasta AA005_Amespa.fasta
```
scp Desktop/data/Luchetti\ data/MITOS_annotated/*.fasta  STUDENTI^pietro.bacconi@137.204.142.152:/DATABIG/pietrobacconi/ncbi_datasets/refseq_mitogenomes/Luchetti_data
```

### downloads

```
mkdir Download_simplifyheaders_MantodeaNCBIdataset/downloads_simplify 

directory((/DATABIG/pietrobacconi/ncbi_datasets/refseq_mitogenomes/Download_MantodeaNCBIdataset)) > cp downloads/*.gene_fasta ../Download_simplifyheaders_MantodeaNCBIdataset/downloads_simplify 


for file in downloads_simplify/*.gene_fasta; do
    base=$(basename "$file" .gene_fasta)
    header=$(head -n 1 "$file" | sed 's/^>//')
    if [[ "$header" == lcl\|NC_* ]]; then
    sed -i -E -e 's/^>lcl\|[^_]+_[0-9]+\.[0-9]+_/>/' -e 's/gene_[0-9]+\_//' -e 's/_\[db_xref=[^]]*\]_//' -e 's/\[location=[^]]*\]_//' -e 's/\[gbkey=[^]]*\]//' "$file"  
    elif [[ "$header" != lcl\|NC_* ]]; then
    sed -i -E 's/^>.*(\[gene=[^]]+\]).*/>\1/' "$file"  
    fi
    sed -E "s/^(>)/\1${base}_/" "$file" > "downloads_simplify/${base}.fasta"
    rm "$file"
done
```

### downloads PCGs

```
mkdir Download_simplifyheaders_MantodeaNCBIdataset/downloads_simplify_PCGs

directory((/DATABIG/pietrobacconi/ncbi_datasets/refseq_mitogenomes/Download_MantodeaNCBIdataset)) > cp downloads_PCGs/*.gene_fasta ../Download_simplifyheaders_MantodeaNCBIdataset/downloads_simplify_PCGs

for file in downloads_simplify_PCGs/*.gene_fasta; do
    base=$(basename "$file" .gene_fasta)
    header=$(head -n 1 "$file" | sed 's/^>//')
    if [[ "$header" == lcl\|NC_* ]]; then
    sed -i -E -e 's/^>lcl\|[^_]+_[0-9]+\.[0-9]+_/>/' -e 's/gene_[0-9]+\_//' -e 's/_\[db_xref=[^]]*\]_//' -e 's/\[location=[^]]*\]_//' -e 's/\[gbkey=[^]]*\]//' "$file"  
    elif [[ "$header" != lcl\|NC_* ]]; then
    sed -i -E 's/^>.*(\[gene=[^]]+\]).*/>\1/' "$file"  
    fi
    sed -E "s/^(>)/\1${base}_/" "$file" > "downloads_simplify_PCGs/${base}.fasta"
    rm "$file"
done
```

### rrns/rrnl

Keeping only the name and the accession

```
mkdir working_directory
cp Mantodea_rrns.fasta Mantodea_rrnl.fasta ../../working_directory/

for file in working_directory/*.fasta; do
sed -i 's/^>\([^:]*\):[^ ]* \([^ ]* [^ ]*\) .*/>\1 \2/' "$file"
sed -i '/^>/ s/ /_/g' "$file"
done 
```

Species names Substitution with code (rrns/rrnl)

```
ln -s /DATABIG/pietrobacconi/ncbi_datasets/refseq_mitogenomes/Download_MantodeaNCBIdataset/downloads_txt/downloads_only_PCGs.txt

awk '
    NR==FNR {
        split($0, a, "_")
        acc = a[1]
        map[acc] = $0
        next
    }
    /^>/ {
        header = substr($0,2)
        split(header, b, "_")
        acc = b[1]
        if (acc in map)
            print ">" map[acc]
        else
            print $0
        next
    }
    { print }
' downloads_only_PCGs.txt Mantodea_rrns.fasta > Mantodea_rrns_renamed.fasta
```
1st error: *Rapttrix_fusca* OM910847.1 (wrong)  Manual final correction headers+sequences *Carrikerella sp.* OM910846.1
2nd error: *Tamolanica_tamolana* (Tamtam header) became *Leptomantella_albella* (Lepalb). Manual correction, pay attention everytime the command is run.   
Final correction adding gene string [gene=rrns]/[gene=rrnl]


```
for file in downloads_simplify_rRNAs/Mantodea_rrnl_renamed.fasta; do
base=$( basename "$file" .fasta)
    sed -i '/^>/ s/$/_\[gene=rrnl\]/' "$file"
done

for file in downloads_simplify_rRNAs/Mantodea_rrns_renamed.fasta; do
base=$( basename "$file" .fasta)
    sed -i '/^>/ s/$/_\[gene=rrns\]/' "$file"
done

mv Mantodea_rrnl_renamed.fasta Mantodea_rrnl_final.fasta
mv Mantodea_rrnl_renamed.fasta Mantodea_rrnl_final.fasta
```


### failed downloads 
```
mkdir Download_MantodeaNCBIdataset/downloads_failed
scp Desktop/data/Failed_downloads/Unverified\ sequences\ MITOS\ annotated/*.fasta  STUDENTI^pietro.bacconi@137.204.142.152:/DATABIG/pietrobacconi/ncbi_datasets/refseq_mitogenomes/Download_MantodeaNCBIdataset/downloads_failed

mkdir Download_simplifyheaders_MantodeaNCBIdataset/downloads_simplify_failed
scp Desktop/data/Failed_downloads/Unverified\ sequences\ MITOS\ annotated/*.fasta  STUDENTI^pietro.bacconi@137.204.142.152:/DATABIG/pietrobacconi/ncbi_datasets/refseq_mitogenomes/Download_simplifyheaders_MantodeaNCBIdataset/downloads_simplify_failed

for file in downloads_simplify_failed/*.fasta; do
    base=$(basename "$file" .fasta)
    sed -i 's/^>.*; *\([^;]*\) */>[gene=\1]/' "$file"
    sed -i "s/^>/>${base}_/" "$file"
    sed -i '/^>/ s/ /_/g' "$file"
done
```

### luchetti data
```
mkdir Download_simplifyheaders_MantodeaNCBIdataset/luchetti_simplify
scp Desktop/data/Luchetti\ data/MITOS_annotated/*.fasta  STUDENTI^pietro.bacconi@137.204.142.152:/DATABIG/pietrobacconi/ncbi_datasets/refseq_mitogenomes/Download_simplifyheaders_MantodeaNCBIdataset/luchetti_simplify

for file in luchetti_simplify/*.fasta; do
    base=$(basename "$file" .fasta)
    sed -i 's/^>.*; *\([^;]*\) */>[gene=\1]/' "$file"
    sed -i "s/^>/>${base}_/" "$file"
    sed -i '/^>/ s/ /_/g' "$file"
done
```

# Grep Headers 20-11-2025

### Grepping mt genes
```
mkdir -p mt_genes_Mantodea
for file in Download_simplifyheaders_MantodeaNCBIdataset/*/*.fasta; do
    awk '
        /^>/ {
            p = 0
            low = tolower($0)

            if (low ~ /gene=(cox1|coxi([^a-z]|$)|coi([^a-z]|$))/)       { out="mt_genes_Mantodea/COX1.fasta"; p=1 }
            else if (low ~ /gene=(cox2|coxii([^a-z]|$)|coii([^a-z]|$))/)       { out="mt_genes_Mantodea/COX2.fasta"; p=1 }
            else if (low ~ /gene=(cox3|coxiii|coiii)/)      { out="mt_genes_Mantodea/COX3.fasta"; p=1 }
            else if (low ~ /gene=(cob|cytb|cyt_b)/)        { out="mt_genes_Mantodea/CYTB.fasta"; p=1 }
            else if (low ~ /gene=(atp6)/)            { out="mt_genes_Mantodea/ATP6.fasta"; p=1 }
            else if (low ~ /gene=(atp8)/)            { out="mt_genes_Mantodea/ATP8.fasta"; p=1 }
            else if (low ~ /gene=(nd1|nad1)/)        { out="mt_genes_Mantodea/ND1.fasta"; p=1 }
            else if (low ~ /gene=(nd2|nad2)/)        { out="mt_genes_Mantodea/ND2.fasta"; p=1 }
            else if (low ~ /gene=(nd3|nad3)/)        { out="mt_genes_Mantodea/ND3.fasta"; p=1 }
            else if (low ~ /gene=(nd4([^a-z]|$)|nad4([^a-z]|$))/)             { out="mt_genes_Mantodea/ND4.fasta"; p=1 }
            else if (low ~ /gene=(nd4l|nad4l)/) { out="mt_genes_Mantodea/ND4L.fasta"; p=1 }
            else if (low ~ /gene=(nd5|nad5)/)        { out="mt_genes_Mantodea/ND5.fasta"; p=1 }
            else if (low ~ /gene=(nd6|nad6)/)        { out="mt_genes_Mantodea/ND6.fasta"; p=1 }
            else if (low ~ /(gene|locus_tag)=(rrn12|rrns|rrnS|s-rrna|.*_mgr01|.*_gr01)/) { out="mt_genes_Mantodea/rrns.fasta"; p=1 }
            else if (low ~ /(gene|locus_tag)=(rrn16|rrnl|rrnL|l-rrna|.*_mgr02|.*_gr02)/) { out="mt_genes_Mantodea/rrnl.fasta"; p=1 }
            else { out="mt_genes_Mantodea/UNMATCHED_genes.fasta"; p=1 }
        }

        p { print >> out }
    ' "$file"
done
```

### Counting tips per gene
```
for file in mt_genes_Mantodea/*.fasta; do

    count=$(grep "^>" "$file" | wc -l)
    echo "Tips number for $file is $count"
done 
```

Until this date 19-11-2025 we have 162 tips from NCBIdataset and 3 Ameles tips, for a total of 165 tips. 
We have 3 missing species in **ND2** fasta file. After review it's clear that *Iris_polystictica*,*Polyspilota_griffinii* and *Otomantis_sp.* have partial genome which there are no ND2 sequences. 


# Dataset Folders Net 11-12-2025 

### Count Families and Subfamilies
```
cut -f 9 Mantodea_NCBIdataset.tsv | sort| uniq | wc -l
cut -f 10 Mantodea_NCBIdataset.tsv | sort | uniq | wc -l
```

In ouput we have **13 Families**, **36 Subfamilies** 

### Counting Genera
```
sed -E 's/_.+$//g' Species_Mantodea.txt > gen.txt 
cat gen.txt | sort | uniq | wc -l      #(94)
mv gen.txt Genera.txt
```

In output we have **94 unique Genera**

## Spltting fasta files into folders
At first we need to copy all the 165 fasta genome files in a single folder which it's called Mantodea_foldersnet. We need to resolve the issue with fasta without rrnl and rrns genes:
```
mv Mantodea_rrnl_final.fasta Mantodea_rrnl_final.fa
mv Mantodea_rrns_final.fasta Mantodea_rrns_final.fa

for file in *.fasta; do
   for rrn_file in *.fa; do
     base=$(basename "$file" .fasta)
     awk -v HEADER=">$base" '
        $0 ~ ("^" HEADER) { p=1; print; next } 
        /^>/{ p=0 } 
        p
     ' "$rrn_file" >> "$file"
   done
done
```
After solving rRNAs issue we start with create the **Folders Network**

```
cut -f 9,10 Mantodea_NCBIdataset.tsv | tail -n +2 - | while IFS=$'\t' read family subfamily; do
  mkdir -p "$family"/"$subfamily"
done
```
Then we **split** all the fasta in their right folders family/subfamily
```
 cut -f 2,9,10 Mantodea_NCBIdataset.tsv | tail -n +2 - | while IFS=$'\t' read code family subfamily; do  
  for file in *.fasta; do
      base=$(basename "$file" .fasta)
      if [[ "$base" == NC_* ]]; then
        name=${base#*_*_}
        elif [[ "$base" != NC_* ]]; then
        name=${base#*_}
      fi
        
      if [[ "$code" == "$name" ]]; then 
      mv "$file" "$family"/"$subfamily"/   
      fi  
  done
done
```
