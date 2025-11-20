# COMANDI GENERICI

Caricare un file da computer al terminale

```
#[Nuovo file modificato caricato sul terminale nella cartella 00_data. Partire dalla directory home del computer]

scp /Desktop/data/Mantodea_dataset.csv STUDENTI^pietro.bacconi@137.204.142.152:/home/STUDENTI/pietro.bacconi/Tirocinio_magistrale/00_data

#[Accedere al server e verificare che il file sia stato caricato. Ora lavorare sulla directory in cui è presente, 00_data.]

cat Mantodea_dataset.csv

#[Caricare un file su DATABIG ha un path diverso. Pertire sempre dalla home del computer]

scp Desktop/data/prova_9AN.txt STUDENTI^pietro.bacconi@137.204.142.152:/DATABIG/pietrobacconi/ncbi_datasets/refseq_mitogenomes

```
Comando per visualizzare un file formato tabella e le sue colonne

```
#[stampa di una colonna]
awk '{ printf $1 }' nome_file.tsv

#[stampa di due colonne separate da tab]
awk '{ printf $1 "\t" $2} nome_file.tsv

```
Comando "grep"
```
#[Visualizzare tutti gli elementi che cominciano con lettera/numero oppure hanno quest parola (esempio con NC)]

grep "NC" AN_Mantodea.txt

#[Contare tutte le righe che hanno quell'elemtento]

grep -c "NC" AN_Mantodea.txt

#[Numerare la riga in cui compare l'elemento]

grep -n "NC" AN_Mantodea.txt

#[Voglio estrarre gli header di un file .fna. Voglio anche contarli]
grep "^>" gene.fna 
grep "^>" gene.fna | wc -l

```
Comando "cut"
```
#[Stampa solo le prime due colonne, usabile come input nei loop]
cut -f1,2 Mantodea_dataset.tsv

```

Scaricare un mitogenoma da NCBI

```
#[entrare nella directory corretta dove si desidera depositare le sequenze]

cd /DATABIG/pietrobacconi/ncbi_datasets/refseq_mitogenoms

#[assicurarsi che l'enviroment di conda sia attivo, sul nostro terminale è stato chiamato downloads]

conda activate downloads
 
#[eseguire il comando del download per il genoma mitocondriale, essendo un mitogenoma useremo il comando 'gene'. Per i nucleari solitamente si usa 'genome'. Il numero di accessp riportato è un esempio, attenzione è importante la dicitura del numero altrimenti il software non lo riconosce correttamente]

datasets download gene accession NC_098567.1

#[Il download restituisce un file zippato che è possibile unzippare e trasferirlo in una cartella a propria scelta (è possibile eseguire due comandi alla volta con l'utilizzo di && o |), il primo comando crea la cartella mentre il secondo unzippa il file e trasferisce il suo contenuto all'interno della cartella creata. Per comodità la cartella viene chiamata con il nome del numero di accesso]

mkdir NC_098567.1/ && unzip NC_098567.1.zip -d NC_098567.1/

```
Modifiche sul dataset ed importo sul terminale. Comandi per contare, ordine alfabetico e stampa

```
#[Modifichiamo la tabella da csv in tsv, togliamo le virgole. Il comando "sed" sostiusce le virgole con tab]

sed 's/,/\t/g' Mantodea_dataset.csv > Mantodea_dataset.tsv
head Mantodea_dataset.tsv

#[Stampiamo solo la prima e la seconda colonna con il comando "awk"]

awk -F'\t' '{print $1, $2}' Mantodea_dataset.tsv

#[Se voglio stampare solo certe righe o escluderne altre devo aggiungere il comando "NR". In questo caso stampa tutte tranne la prima, che è l'header]

awk -F'\t' 'NR > 1 {print $1, $2}' Mantodea_dataset.tsv

#[Se voglio mettere le righe in ordine alfabetico/numerico devo usare il comando "sort --unique". La | mi permette di svolgere entrambi i comandi]

awk -F'\t' 'NR > 1 {print $1, $2}' Mantodea_dataset.tsv | sort --unique

#[Se voglio contare quanti elementi ci sono per colonna devo usare il comando "wordcount"]

awk -F'\t' '{print $1, $2}' Mantodea_dataset.tsv | sort --unique | wc -l

```
SACRICARE UN DATASET (Accession Numbers + Species Name)     Primo tentativo 22/10/2025

```
conda activate downloads

paste prova_9AN.txt prova_9species.txt | while IFS=$'\t' read acc species; 
do
    echo "Downloading $acc for $species..."

    safe_species=$(echo "$species" | tr ' ' '_' | tr -d '()')

    datasets download gene accession "$acc" --filename "${acc}.zip"

    unzip -o "${acc}.zip" -d "$safe_species"

    echo "✓ Downloaded and extracted to $safe_species"
done

```
- Il risultato di paste è preso in input dal while tramite la pipe.
- Comando paste unisce i due file in modo da poter leggere contemporaneamente l'AN e il nome della specie. La cosa FONDAMENTALE è assicurarsi che i due txt abbiano lo stesso ordine.
- Nella prima riga dopo il 'do' sto salvando nella nuova variabile safe_species il nome da dare in output tramite echo preso dalla variabile $species, nella quale viene copiato al suo interno il nome di una specie dalla lista del file .txt per ogni iterazione (grazie al comando read). La seconda parte di comando è essenziale per correggere i nomi e togliere potenziali caratteri problematici.
- Se non esiste la directory non è necessario usare mkdir, il comando unzip -d crea la directory se non esiste.

SCRIPT DOWNLOAD DI UN DATASET MIGLIORATO

```
conda activate downloads

cut -f1,2 Mantodea_dataset.tsv |tail -n +2 Mantodea_dataset.tsv | while IFS=$'\t' read acc species; 
do

echo "Downloading $acc for $species..."

folder_name="${species}_${acc}"

datasets download gene accession "$acc" --include gene,rna,cds,protein,5p-utr,3p-utr,product-report --filename "${species}_${acc}.zip"

unzip -o "${species}_${acc}.zip" -d "$folder_name"

mv "$folder_name"/ncbi_dataset/data/gene.fna "$folder_name"/ncbi_dataset/data/"${species}_${acc}_genes.fna"
mv "$folder_name"/ncbi_dataset/data/protein.faa "$folder_name"/ncbi_dataset/data/"${species}_${acc}_protein.faa"
mv "$folder_name"/ncbi_dataset/data/data_report.jsonl "$folder_name"/ncbi_dataset/data/"${species}_${acc}_data_report.jsonl"
mv "$folder_name"/ncbi_dataset/data/dataset_catalog.jsonl "$folder_name"/ncbi_dataset/data/"${species}_${acc}_dataset_catalog.jsonl"
mv "$folder_name"/ncbi_dataset/data/product_report.jsonl "$folder_name"/ncbi_dataset/data/"${species}_${acc}_product_report.jsonl"

find "$folder_name" -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \) | while read fasta; do
        sed -i '/^>/ s/ /_/g' "$fasta"
    done

mv "${species}_${acc}.zip" Mantodea_zip/

echo "✓ Downloaded and extracted to $folder_name"

done

```

-Aggiunto il comando cut per considerare in input solo le prime due colonne
-Aggiunto --include per scaricare correttamente tutto il contenuto del mitogenoma
-Aggiunto 5 line per rinominare tutti file ottenuti dall'unzip 
-Aggiunto comando eliminazione degli spazi negli header

DOWNLOAD NOTES 4/11/2025 
```
mkdir -p downloads_only_13 correct_downloads
set -uo pipefail

while IFS=$'\t' read -r acc;

 do
    [ -z "$acc" ] && continue
    echo "Processing ${acc}.."
    if ! esearch -db nucleotide -query "${acc}" </dev/null | efetch -format gene_fasta > "${acc}.gene_fasta"; then
       echo "$acc" >> failed_downloads.txt
       echo "Download "$acc" failed"
       continue
    fi
    [ -s "${acc}.gene_fasta" ] || {
    echo "${acc}.gene_fasta is empty"
    echo "$acc" >> failed_downloads.txt
    rm "${acc}.gene_fasta"
    continue
    } 
       count=$( (grep -c "^>" "${acc}.gene_fasta" 2>/dev/null) || echo 0 )
    
    if  [ "$count" -le 13 ]; then
         
         sed -i '/^>/ s/ /_/g' ${acc}.gene_fasta
         echo "$acc" >> downloads_only_13.txt
         echo "$acc has only $count sequences"
         mv "${acc}.gene_fasta" downloads_only_13/

       else
       
       sed -i '/^>/ s/ /_/g' ${acc}.gene_fasta
       echo "$acc downloaded correctly with $count sequences"
       echo "$acc" >> correct_downloads.txt
       mv "${acc}.gene_fasta" correct_downloads/

    fi

done < AN_Mantodea.txt

```
#DOWNLOADS FALLITI (7):
OZ336339.1  assembly   Bolivaria brachyptera
PP438765.1  unverified Ceratocrania macra
PP438767.1  unverified Chlidonoptera lestoni
MN267041.1  unverified Mekongomantis quinquespinosa
PP438763.1  unverified Astyliasula basinigra
PP438766.1  unverified Ceratomantis yunnanensis
PP438770.1  unverified Panurgica fratercula 

#DOWNLOADS SOLO CDS 12_13 (79)

#DOWNLOAD SUCCESFULLY (76)



DOWNLOAD NOTES 11/11/2025
```
mkdir -p downloads_only_13 correct_downloads
set -uo pipefail
 
while IFS=$'\t' read species acc;

 do
    [ -z "$acc" ] && continue
    echo "Processing ${acc} for ${species}..."
    if ! esearch -db nucleotide -query "${acc}" </dev/null | efetch -format gene_fasta > "${acc}_${species}.gene_fasta"; then
       echo "${acc}_${species}" >> failed_downloads.txt
       echo "Download "${acc}_${species}" failed"
       continue
    fi
    [ -s "${acc}_${species}.gene_fasta" ] || {
    echo "${acc}_${species}.gene_fasta is empty"
    echo "${acc}_${species}" >> failed_downloads.txt
    rm "${acc}_${species}.gene_fasta"
    continue
    } 
       count=$( (grep -c "^>" "${acc}_${species}.gene_fasta" 2>/dev/null) || echo 0 )
    
    if  [ "$count" -le 13 ]; then
         
         sed -i '/^>/ s/ /_/g' ${acc}_${species}.gene_fasta
         echo "${acc}_${species}" >> downloads_only_13.txt
         echo "${acc}_${species} has only $count sequences"
         mv "${acc}_${species}.gene_fasta" downloads_only_13/

       else
       
       sed -i '/^>/ s/ /_/g' ${acc}_${species}.gene_fasta
       echo "${acc}_${species} downloaded correctly with $count sequences"
       echo "${acc}_${species}" >> correct_downloads.txt
       mv "${acc}_${species}.gene_fasta" correct_downloads/
    fi
done < C1_C2_Mantodea_dataset.tsv 
```



volendo potrebbe funzionare:

```
<(cut -f 1,2 Mantodea_dataset.tsv) | tail -n +2 Mantodea_dataset.tsv | while IFS=$'\t' read species acc;
do
......
......
......
done

```
WRITING NEW HEADERS 12/11/2025

```
#downloads
mkdir downloads_simplify
for file in downloads/*.gene_fasta; do
    base=$(basename "$file" .gene_fasta)
    sed -E -e 's/^>lcl\|[^_]+_[0-9]+\.[0-9]+_/>/' -e 's/gene_[0-9]+\_//' -e 's/_\[db_xref=[^]]*\]_//' -e 's/\[location=[^]]*\]_//' -e 's/\[gbkey=[^]]*\]//' "$file" > "downloads_simplify/${base}.fasta"  
    sed -i -E "s/^(>)/\1${base}_/" "downloads_simplify/${base}.fasta"
done


#downloads_PCGs_13
mkdir downloads_PCGs_13_simplify
for file in downloads_PCGs_13/*.gene_fasta; do
    base=$(basename "$file" .gene_fasta)
    sed -E -e 's/^>lcl\|[^_]+_[0-9]+\.[0-9]+_/>/' -e 's/gene_[0-9]+\_//' -e 's/\[location=[^]]*\]_//' -e 's/\[gbkey=[^]]*\]//' "$file" > "downloads_PCGs_13_simplify/${base}.fasta"  
    sed -i -E "s/^(>)/\1cd ..${base}_/" "downloads_PCGs_13_simplify/${base}.fasta"
done



#downloads_rRNAs
mkdir downloads_rRNAs_simplify
cp downloads_rRNAs/Mantodea_16S.fasta downloads_rRNAs/Mantodea_12S.fasta downloads_rRNAs_simplify/
for file in downloads_rRNAs_simplify/*.fasta; do
    base=$(basename "$file" .fasta)
    sed -i 's/^>\([^:]*\):[^ ]* \([^ ]* [^ ]*\) .*/>\1 \2/' "$file" 
    sed -i '/^>/ s/ /_/g' "$file"
    if [ "$base" = "Mantodea_16S" ]; then
    sed -i '/^>/ s/$/_\[gene=16S\]/' "$file"
    else 
    sed -i '/^>/ s/$/_\[gene=12S\]/' "$file"
    fi
done


#14/11/2025
#downloads failed 

mkdir downloads_failed_simplify
scp Desktop/data/Unverified\ sequences\ MITOS\ annotated/* STUDENTI^pietro.bacconi@137.204.142.152:/DATABIG/pietrobacconi/ncbi_datasets/refseq_mitogenomes/Downloads_Mantodea_simplify/downloads_failed_simplify 

for file in downloads_failed_simplify/*.fasta; do
    base=$(basename "$file" .fasta)
    sed -i 's/^>.*; *\([^;]*\) */>[gene=\1]/' "$file"
    sed -i "s/^>/>${base}_/" "$file"
    sed -i '/^>/ s/ /_/g' "$file"
done

#14/11/2025
#luchetti_data_3_species
mkdir luchetti_data_simplify
scp Desktop/data/Luchetti\ data/MITOS_annotated/*  STUDENTI^pietro.bacconi@137.204.142.152:/DATABIG/pietrobacconi/ncbi_datasets/refseq_mitogenomes/Downloads_Mantodea_simplify/luchetti_data_simplify/

for file in luchetti_data_simplify/*.fasta; do
    base=$(basename "$file" .fasta)
    sed -i 's/^>.*; *\([^;]*\) */>[gene=\1]/' "$file"
    sed -i "s/^>/>${base}_/" "$file"
    sed -i '/^>/ s/ /_/g' "$file"
done
```

CREATING GENES FILES (GREPPING HEADERS+SEQUENCE) 14-11-2025

```
mkdir -p mt_genes_Mantodea
for file in Downloads_Mantodea_simplify/*/*.fasta; do
    awk '
        /^>/ {
            p = 0
            low = tolower($0)

            if (low ~ /gene=(cox1|coxi([^a-z]|$)|coi([^a-z]|$))/)       { out="mt_genes_Mantodea/COX1_genes.fasta"; p=1 }
            else if (low ~ /gene=(cox2|coxii([^a-z]|$)|coii([^a-z]|$))/)       { out="mt_genes_Mantodea/COX2_genes.fasta"; p=1 }
            else if (low ~ /gene=(cox3|coxiii|coiii)/)      { out="mt_genes_Mantodea/COX3_genes.fasta"; p=1 }
            else if (low ~ /gene=(cob|cytb|cyt_b)/)        { out="mt_genes_Mantodea/CYTB_genes.fasta"; p=1 }
            else if (low ~ /gene=(atp6)/)            { out="mt_genes_Mantodea/ATP6_genes.fasta"; p=1 }
            else if (low ~ /gene=(atp8)/)            { out="mt_genes_Mantodea/ATP8_genes.fasta"; p=1 }
            else if (low ~ /gene=(nd1|nad1)/)        { out="mt_genes_Mantodea/ND1_genes.fasta"; p=1 }
            else if (low ~ /gene=(nd2|nad2)/)        { out="mt_genes_Mantodea/ND2_genes.fasta"; p=1 }
            else if (low ~ /gene=(nd3|nad3)/)        { out="mt_genes_Mantodea/ND3_genes.fasta"; p=1 }
            else if (low ~ /gene=(nd4([^a-z]|$)|nad4([^a-z]|$))/)             { out="mt_genes_Mantodea/ND4_genes.fasta"; p=1 }
            else if (low ~ /gene=(nd4l|nad4l)/) { out="mt_genes_Mantodea/ND4L_genes.fasta"; p=1 }
            else if (low ~ /gene=(nd5|nad5)/)        { out="mt_genes_Mantodea/ND5_genes.fasta"; p=1 }
            else if (low ~ /gene=(nd6|nad6)/)        { out="mt_genes_Mantodea/ND6_genes.fasta"; p=1 }
            else if (low ~ /(gene|locus_tag)=(rrn12|rrns|12s|.*_mgr01)/) { out="mt_genes_Mantodea/rrns_genes.fasta"; p=1 }
            else if (low ~ /(gene|locus_tag)=(rrn16|rrnl|16s|.*_mgr02)/) { out="mt_genes_Mantodea/rrnl_genes.fasta"; p=1 }
            else { out="mt_genes_Mantodea/UNMATCHED_genes.fasta"; p=1 }
        }

        p { print >> out }
    ' "$file"
done



#COUNTING TIPS PER GENE
for file in mt_genes_Mantodea/*.fasta; do

count=$(grep "^>" "$file" | wc -l)
echo "Tips number for $file is $count"

done 
```
OUTPUT:
Tips number for mt_genes_Mantodea/ATP6_genes.fasta is 165
Tips number for mt_genes_Mantodea/ATP8_genes.fasta is 165
Tips number for mt_genes_Mantodea/COX1_genes.fasta is 165
Tips number for mt_genes_Mantodea/COX2_genes.fasta is 165
Tips number for mt_genes_Mantodea/COX3_genes.fasta is 165
Tips number for mt_genes_Mantodea/CYTB_genes.fasta is 165
Tips number for mt_genes_Mantodea/ND1_genes.fasta is 165
Tips number for mt_genes_Mantodea/ND2_genes.fasta is 160
Tips number for mt_genes_Mantodea/ND3_genes.fasta is 165
Tips number for mt_genes_Mantodea/ND4_genes.fasta is 165
Tips number for mt_genes_Mantodea/ND4L_genes.fasta is 165
Tips number for mt_genes_Mantodea/ND5_genes.fasta is 165
Tips number for mt_genes_Mantodea/ND6_genes.fasta is 165
Tips number for mt_genes_Mantodea/rrnl_genes.fasta is 156
Tips number for mt_genes_Mantodea/rrns_genes.fasta is 156
Tips number for mt_genes_Mantodea/UNMATCHED_genes.fasta is 1917


## FINDING MISSING SPECIES in grep process  17-11-2025 / 18-11-2025
 
#ND2
```
mkdir missing_headers_ND2
cp ND2_genes.fasta missing_headers_ND2/

#SERVE PER CORREGGERE GLI HEADER ED ELIMINARE LE LCL DI TROPPO 
sed -i -E 's/lcl[^.]*\.1_//' file.fasta

grep "^>" ND2_genes.fasta > ND2_output.fasta  
sed -i -E 's/^>.*_([A-Za-z]+_[A-Za-z]+).*/\1/' ND2_output.fasta
sort -u ND2_output.fasta > species_in_ND2_fasta.txt 
sort -u ../../Species_Mantodea.txt > species_list_sorted.txt
diff species_in_ND2_fasta.txt species_list_sorted.txt > missing.txt
```

Output missing.txt per ND2:
Tenodera_aridifolia_brevicollis 
Decimiana_sp. //check su NCBI, ND2 PRESENTE su GenBank ma NON SCARICATA, BISOGNA AGGIUNGERLE A MANO
Deroplatys_desiccata_2                
Hierodula_sp._2q
Humbertiella_nada_2
Iris_polystictica //chack su NCBI, ND2 assente partial genomecd
Mesopteryx_alata_2
Phyllothelys_breve_2
Polyspilota_griffinii //check su NCBI, ND2 assente partial genome
Popa_spurca_spurca
Pseudempusa_pinnapavonis_2
Spilomantis_occipitalis_2
Statilia_sp._2
Arria_pallida_2
Taumantis_sigiana // check su NCBI, ND2 PRESENTE su GenBank ma NON SCARICATA, BISOGNA AGGIUNGERLE A MANO
Otomantis_sp. //chack su NCBI, ND2 assente partial genome


#rrns
``` 
mkdir missing_headers_rrns
cp rrns_genes.fasta missing_headers_rrns/

sed -i -E 's/lcl[^.]*\.1_//' file.fasta 

grep "^>" rrns_genes.fasta > rrns_output.fasta  
sed -i -E 's/^>.*_([A-Za-z]+_[A-Za-z]+).*/\1/' rrns_output.fasta
sort -u rrns_output.fasta > species_in_rrns_fasta.txt 
sort -u ../../Species_Mantodea.txt > species_list_sorted.txt
diff species_in_rrns_fasta.txt species_list_sorted.txt > missing.txt
```
utput:

Carrikerella_sp.  //check su NCBI rrns PRESENTE, trovato errore nell'accession ORA CORRETTO 

Gonypeta_sp.  //check su NCBI rrns PRESENTE, pattern di ricerca diverso, "s-rRNA" (l-rRNA per 16S)


Popa_spurca_spurca //check c'é, il seconda spurca è venuto tagliato 

Stenotoxodera_porioni //check su NCBI rrns PRESENTE, HEADER SBAGLIATO nei download, CORRETTI A MANO, RIFARE grep
Theopropus_elegans  //check su NCBI rrns PRESENTE, HEADER SBAGLIATO nei download, CORRETTI A MANO, RIFARE grep


Anaxarcha_zhengi       //check su NCBI rrns PRESENTE, pattern di ricerca diverso "gr01" (gr02 per 16S)
Creobroter_gemmatus    //check su NCBI rrns PRESENTE, pattern di ricerca diverso "gr01" (gr02 per 16S)
Hierodula_formosana    //check su NCBI rrns PRESENTE, pattern di ricerca diverso "gr01" (gr02 per 16S)
Mantis_religiosa       //check su NCBI rrns PRESENTE, pattern di ricerca diverso "gr01" (gr02 per 16S)
Tenodera_sinensis      //check su NCBI rrns PRESENTE, pattern di ricerca diverso "gr01" (gr02 per 16S)


## FINDING MISSING SPECIES in grep process 2°run 20-11-2025

File .txt with code names
```
(~/Tirocinio_magistrale/00_data/Mantodea_data) -> cut -f2 Mantodea_NCBIdataset.tsv > Code_MantodeaNCBIdataset.txt
(/DATABIG/pietrobacconi/ncbi_datasets/refseq_mitogenomes/mt_genes_Mantodea) -> ln -s /home/STUDENTI/pietro.bacconi/Tirocinio_magistrale/00_data/Mantodea_data/Code_MantodeaNCBIdataset.txt
```

Manual: Delete header and had Ameand Ameser Amespa
```
nano Code_MantodeaNCBIdataset.txt
```

Finding missing species
```
grep "^>" rrns.fasta > rrns_output.fasta  
sed -E -i 's/^>NC[^_]*_[^_]*_([^_]*).*/>\1/; s/^>[^_]*_([^_]*).*/>\1/' rrns_output.fasta
sort rrns_output.fasta > species_in_rrns_fasta.txt 
sort Code_MantodeaNCBIdataset.txt > species_list_sorted.txt
diff species_in_rrns_fasta.txt species_list_sorted.txt > missing.txt

grep "A" missing.txt | less
grep "B" missing.txt | less
grep "C" missing.txt | less
...
```

### Errors:

- Gonypeta_sp.: (downloads_simplify) grep didn't work 
- Hierodula_sp.: (downloads_simplify) grep didn't work 
- Staspe2: (downloads_simplify) grep didn't work -> pattern s-rRNA -> if condition corretta in s-rrna
- Leptomantella albella è venuto doppio???? Find error: don't know why Tamtam became Lepalb during the simplify process CORRECT VERSION -> NC_007702.1_Tamtam_[gene=rrnl]

Run once again the loop with pattern correction, right result

