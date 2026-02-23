#!/bin/bash

# Mapping short reads on assembly. Mutate SAM into BAM
minimap2 -ax sr --MD -t 6 node_18_runall.fasta Amedum_1_paired.fastq  Amedum_2_paired.fastq > Amedum_raw_sr.sam
samtools view -Sb Amedum_raw_sr.sam > Amedum_raw_sr.bam
rm Amedum_raw_sr.sam
samtools sort -@6 -o Amedum_raw_sr_sorted.bam Amedum_raw_sr.bam
samtools index Amedum_raw_sr_sorted.bam
rm Amedum_raw_sr.bam
