COMANDI GENERICI

Comando per visualizzare un file formato tabella e le sue colonne


```
#[stampa di una colonna]
awk '{ printf $1 }' nome_file.tsv

#[stampa di due colonne separate da tab]
awk '{ printf $1 "\t" $2} nome_file.tsv

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
#[Nuovo file modificato caricato sul terminale nella cartella 00_data. Partire dalla directory home del computer]

scp /Desktop/data/Mantodea_dataset.csv STUDENTI^pietro.bacconi@137.204.142.152:/home/STUDENTI/pietro.bacconi/Tirocinio_magistrale/00_data

#[Accedere al server e verificare che il file sia stato caricato. Ora lavorare sulla directory in cui è presente, 00_data.]

cat Mantodea_dataset.csv

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
