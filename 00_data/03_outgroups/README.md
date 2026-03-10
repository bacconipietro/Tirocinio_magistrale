# Selected Outgroups

| Species | Code | Accession | Family |
|---|---|---|---|
| Reticulitermes_santonensis | Retsan | NC_009499.1 | Termitoidae |
| Mastotermes_darwiniensis | Masdar | NC_018120.1 | Termitoidae |
| Coptotermes_formosanus | Copfor | NC_015800.1 | Termitoidae |
| Periplaneta_brunnea | Perbru | NC_039940.1 | Blattidae |
| Blatella_germanica | Blager | NC_012901.1 | Ectobiidae |
| Diploptera_punctata | Dippun | NC_082459.1 | Blaberidae |
| Eupolyphaga_sinensis | Eupsin | NC_014274.1 | Corydiidae |
| Cryptocercus_meridianus | Crymer | NC_037496.1 | Cryptocercidae |

## Download Outgroups 10-03-2026
```bash
mkdir -p downloads downloads_PCGs downloads_failed downloads_txt
cut -f2,3 selected_outgroups.tsv | while IFS=$'\t' read species acc;
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
downloads_PCGs, Manual adding and modifying rRNAs:
+ NC_015800.1_Copfor.gene_fasta
+ NC_009499.1_Retsan.gene_fasta  
+ NC_012901.1_Blager.gene_fasta  
+ NC_014274.1_Eupsin.gene_fasta  

#### Flip reversed sequences
```bash
bash reverse_downloads_complement.sh /DATABIG/pietrobacconi/ncbi_datasets/refseq_mitogenomes/00_data/00_NCBI/00_Mitochondrial/Download_OutgroupNCBIdataset/00_download/downloads -v -o /DATABIG/pietrobacconi/ncbi_datasets/refseq_mitogenomes/00_data/00_NCBI/00_Mitochondrial/Download_OutgroupNCBIdataset/01_flipped/flipped_downloads

bash reverse_downloads_complement.sh /DATABIG/pietrobacconi/ncbi_datasets/refseq_mitogenomes/00_data/00_NCBI/00_Mitochondrial/Download_OutgroupNCBIdataset/00_download/downloads_PCGs -v -o /DATABIG/pietrobacconi/ncbi_datasets/refseq_mitogenomes/00_data/00_NCBI/00_Mitochondrial/Download_OutgroupNCBIdataset/01_flipped/flipped_downloads_PCGs
```
#### Simplify headers
```bash
for file in *.gene_fasta; do
    base=$(basename "$file" .gene_fasta)
    header=$(head -n 1 "$file" | sed 's/^>//')
    if [[ "$header" == lcl\|NC_* ]]; then
    sed -i -E -e 's/^>lcl\|[^_]+_[0-9]+\.[0-9]+_/>/' -e 's/gene_[0-9]+\_//' -e 's/_\[db_xref=[^]]*\]_//' -e 's/\[location=[^]]*\]_//' -e 's/\[gbkey=[^]]*\]//' "$file"  
    elif [[ "$header" != lcl\|NC_* ]]; then
    sed -i -E 's/^>.*(\[gene=[^]]+\]).*/>\1/' "$file"  
    fi
    sed -E "s/^(>)/\1${base}_/" "$file" > "${base}.fasta"
    rm "$file"
done
```
