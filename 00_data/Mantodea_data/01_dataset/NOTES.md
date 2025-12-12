# Dataset notes

Starting row numbers 163, with header 
(wh = without header)

```
#Species number without sorting wh == 162
awk -F'\t' 'NR > 1 {print $1}' Mantodea_dataset.tsv | wc -l

#Unique species number sorted wh == 153. There are 9 double species

awk -F'\t' 'NR > 1 {print $1}' Mantodea_dataset.tsv | sort --unique | wc -l

#Accessions number without sorting wh == 162]
awk -F'\t' 'NR > 1 {print $2}' Mantodea_dataset.tsv |  wc -l

#Unique accessions number without sorting wh == 162]
awk -F'\t' 'NR > 1 {print $2}' Mantodea_dataset.tsv | sort --unique | wc -l

```
In this dataset there are 162 unique AN, 153 of these are unique species meanwhile 9 are double.



AGG. 14/11/2025
We add three new species (AA001,AA003,AA005), reaching 165 total dataset examples.  


