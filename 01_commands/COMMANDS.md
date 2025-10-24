COMANDI GENERICI

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

