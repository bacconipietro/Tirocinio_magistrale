#Comando_awk

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
 
#[eseguire il comando del download per il genoma mitocondriale, essendo un mitogenoma useremo il comando 'gene'.
Per i nucleari solitamente si usa 'genome'.
Il numero di accessp riportato è un esempio,
attenzione è importante la dicitura del numero altrimenti il software non lo riconosce correttamente]

datasets download gene accession NC_098567.1

#[Il download restituisce un file zippato che è possibile unzippare e trasferirlo in una cartella a propria scelta (è possibile eseguire due comandi alla volta con l'utilizzo di && o |), il primo comando crea la cartella mentre il secondo unzippa il file e trasferisce il suo contenuto all'interno della cartella creata. Per comodità la cartella viene chiamata con il nome del numero di accesso]

mkdir NC_098567.1/ && unzip NC_098567.1.zip -d NC_098567.1/

```
