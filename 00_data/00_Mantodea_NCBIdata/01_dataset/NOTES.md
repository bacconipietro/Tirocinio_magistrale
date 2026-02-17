# Dataset notes

Starting row numbers 163, with header 
(wh = without header)

```bash
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

# 09-02-2026
We add Paratoxodera gigliotosi MG888457.1, reaching **163** total NCBI samples 

# 17-02-2026
We add 7 new species, reaching **170** total NCBI samples
Ambivia_parapopa	            MG888436.1
Amorphoscelis_sp.	            MG888437.1
Hymenopus_coronatoides	      MG888449.1
Leptomantella_lactea	        MG888453.1
Metallyticus_splendidus	      MG888455.1
Parablepharis_kuhlii_asiatica	MG888456.1
Sinomiopteryx_guangxiensis	  MG888464.1
