#Notes sul dataset

Numero di righe iniziale 163, header compreso 
(wh = wuthout header)

```
#Numero specie no ordine wh == 162]
awk -F'\t' 'NR > 1 {print $1}' Mantodea_dataset.tsv | wc -l

#Numero specie ordine unvico wh == 153. Ci sono 9 doppioni]

awk -F'\t' 'NR > 1 {print $1}' Mantodea_dataset.tsv | sort --unique | wc -l

#Numero AN no ordine wh == 162]
awk -F'\t' 'NR > 1 {print $2}' Mantodea_dataset.tsv |  wc -l

#Numero AN ordine univoco wh == 162]
awk -F'\t' 'NR > 1 {print $2}' Mantodea_dataset.tsv | sort --unique | wc -l

```
Nel dataset sono presenti 162 Acession Numbers univoci, di questi 153 sono di Specie univoche mentre 9 sono doppioni



AGG. in data 14/11/2025

Nel dataset sono presenti 162 Acession Numbers univoci, di questi 153 sono di Specie univoche mentre 9 sono doppioni.
A questi si aggiungono altri 3 campioni (AA001,AA003,AA005) per un totale di 165 tips. 


