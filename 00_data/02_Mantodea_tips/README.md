## Count Families and Subfamilies 11-12-2025  
```
cut -f 9 Mantodea_NCBIdataset.tsv | sort| uniq | wc -l
cut -f 10 Mantodea_NCBIdataset.tsv | sort | uniq | wc -l
```

In ouput we have **13 Families**, **36 Subfamilies** 

### Counting Genera
```
sed -E 's/_.+$//g' Species_Mantodea.txt > gen.txt 
cat gen.txt | sort | uniq | wc -l      #(94)
mv gen.txt Genera.txt
```

In output we have **94 unique Genera**

-----

## Count update 04-03-2026

```bash
#Sorting taxonomy file:
tail -n +2 Taxonomy.tsv | sort -k1,1 -k2,2 | uniq | cat <(head -n 1 Taxonomy.tsv) - > taxonomy_sorted.tsv

#counting unique families from the same file
tail -n +2 taxonomy_sorted.tsv | cut -f1 | sort -u | wc -l
```

Update: Families **19**, Subfamilies **33**, Genera **98**. 
