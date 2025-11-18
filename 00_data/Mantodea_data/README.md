##Download NCBI dataset 18/11/2025

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


