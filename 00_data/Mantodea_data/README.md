# Download NCBI dataset 18/11/2025

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
One error: Rapttrix_fusca OM910847.1 X -> Carrikerella sp. OM910846.1 //Manual final correction headers+sequences

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




