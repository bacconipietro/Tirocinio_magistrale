#Comando_awk

Comando per visualizzare un file formato tabella e le sue colonne


```
#[stampa di una colonna]
awk '{ printf $1 }' nome_file.tsv

#[stampa di due colonne separate da tab]
awk '{ printf $1 "\t" $2} nome_file.tsv

```


