# Trimming short Reads 16-12-2025

#### Installing trimmomatic
```bash
wget https://github.com/usadellab/Trimmomatic/releases/download/v0.40/Trimmomatic-0.40.zip
unzip Trimmomatic-0.40.zip
java --version #check java version
conda create --name assembly
conda activate assembly
conda install -c bioconda trimmomatic
```

#### Run trimmomatic

```bash
#general 
trimmomatic PE -threads 20 -phred33 ../00_data/<SPECIE_READS_CODE>_1.fastq.gz ../00_data/<SPECIE_READS_CODE>_2.fastq.gz <SPECIES_CODE>_1_paired.fastq <SPECIES_CODE>_1_unpaired.fastq <SPECIES_CODE>_2_paired.fastq <SPECIES_CODE>_2_unpaired.fastq ILLUMINACLIP:../trimmomatic-0.40/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> 00_stats/<SPECIES_CODE>_stats_trimmomatic.log 

#example of Pseudoyersinia betancuriae
trimmomatic PE -threads 20 -phred33 ../00_data/16-CI1f_1.fastq.gz ../00_data/16-CI1f_2.fastq.gz Psebet_1_paired.fastq Psebet_1_unpaired.fastq Psebet_2_paired.fastq Psebet_2_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> Psebet_stats_trimmomatic.log 
```

-----

# Assembly with SPAdes and lcWGS
 
## SPAdes 17-12-2025

#### Install and test 
```bash
conda activate assembly
wget https://github.com/ablab/spades/releases/download/v4.2.0/SPAdes-4.2.0-Linux.tar.gz
conda install -c bioconda SPAdes

#test run
spades.py --test
```

#### Subsamples short reads

We subsample reads data to 2Mbp file sequence. If there's no positive result after assembling it's possible to rise subsample content (e.g. 4Mbp, 8Mbp, 10Mbp, 12Mbp).
```bash

#general
for file in 01_trimmed/*/*_<FORWARD OR REVERSE>_paired.fastq; 
do  
base=$(basename "$file" .fastq) 
seqtk sample -s100 "$file" <READS_NUMBER> > "${base}_<FORWARD OR REVERSE>_<READS_NUMBER>Mbp.fq" 
done


#Subsample of 2Mbp
for file in 01_trimmed/*/*_1_paired.fastq; 
do  
base=$(basename "$file" .fastq) 
seqtk sample -s100 "$file" 2000000 > "${base}_1_2Mbp.fq" 
done


for file in 01_trimmed/*/*_2_paired.fastq; 
do  
base=$(basename "$file" .fastq) 
seqtk sample -s100 "$file" 2000000 > "${base}_2_2Mbp_.fq" 
done


#counting reads number to check if the for worked
zgrep -c "^@" <SPECIES_CODE>_1_2Mbp.fq 
```

#### Run SPAdes

Assembly output show `contigs.fasta`. After check contig length (15-17kbp) and coverage (>=10x) extract candidate contigs and BLAST them to see if are what you are looking for.   
```bash
#general
spades.py --only-assembler -t 8 -1 <SPECIES_CODE>_1_2Mbp.fq -2 <SPECIES_CODE>_2_2Mbp.fq -o <SPECIES_CODE>

#extraction candidate contigs
awk '/^>NODE_<NODE_NUMBER>_/{print;p=1;next} /^>/{p=0} p' contigs.fasta > node_<NODE_NUMBER>.fasta    
```
After blasting, going on [MITOS2](https://usegalaxy.org/root?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fiuc%2Fmitos2%2Fmitos2%2F2.1.3%20galaxy0) to annotate mitogenomes, avaible on __GALAXY__, open source web-based platform for data intensive biomedical research.

MITOS2 options:
+ Invertebrate (5)
+ RefSeq63 Metazoa
+ output:nucleotide FASTA
+ Advances options: Final overlap 50 - fragment overlap 10 - flag annotate best options
+ Advances options for protein: select first two flags down below

#### Flip possible reversed sequences:
```bash
bash /home/STUDENTI/pietro.bacconi/Tirocinio_magistrale/99_scripts/reverse_downloads_failed.sh *.fasta
```

#### Simplify headers
```bash
for file in work/*.fasta; do
    base=$(basename "$file" .fasta)
    sed -i 's/^>.*; *\([^;]*\) */>[gene=\1]/' "$file"
    sed -i "s/^>/>${base}_/" "$file"
    sed -i '/^>/ s/ /_/g' "$file"
done
```
-----

## lcWGS 30-01-2026

#### Theory:

Your script is designed to run a "parameter sweep"—it tries to assemble the genome using increasing numbers of reads (50k, 100k, etc.) to find the minimum required for a good assembly.
It's data-driven: The assembler (MitoZ) looks for an overlap between the start and end of the sequence. If the coverage drops at the ends or the sequence is repetitive, it cannot confidently join them, so it outputs a linear sequence.
The Annotation consequence: Because the output file is technically linear, the Annotation step (which gives you the warning) treats it as a straight line. This is why your tRNAs are missing. If a gene sits exactly across the "break" point of the circle, the linear annotator cannot see it.

+ Run `run_mt_sweep_single.sh`:
```bash
ln -s /home/STUDENTI/pietro.bacconi/Tirocinio_magistrale/99_scripts/lcWGS/lcWGS_env.from_history.yaml
ln -s /home/STUDENTI/pietro.bacconi/Tirocinio_magistrale/99_scripts/lcWGS/run_mt_sweep_single.sh

mamba create -f lcWGS_env.from_history.yaml -n lcWGS
mamba activate lcWGS

bash run_mt_sweep_single.sh -1 16-CI1f_1.fastq.gz  -2 16-CI1f_2.fastq.gz -s Pseben -r sweep_single_pseben
```

Once canditate contigs are reached, it's necessary to merge Annotation information with respective contig. 

```bash 
#cp cordinates in > positions_Psebet.txt

less mt_summary.txt

#redoo for all lcWGS outputs
awk '{print $1"\t"($2-1)"\t"$3"\t"$7"\t.\t"$5}' positions_Psebet.txt > positions_Psebet.bed
seqkit subseq --bed positions_Psebet.bed 16CI1f_Psebet.fasta -o annotated_16CI1f_Psebet.fasta
```

#### Flip possible reversed seuquences
```bash
bash /home/STUDENTI/pietro.bacconi/Tirocinio_magistrale/99_scripts/reverse_lcWGS.sh *.fasta 
```

#### Simplify headers
```bash
for file in *.fasta; do
    base=$(basename "$file" .fasta)
    sed -i -E "s/>NODE_[0-9]+_[0-9]+-[0-9]+:[+-] (.*)/>${base}_[gene=\1]/" "$file"
done
```

## Resolving Amedum 17-02-2026
Every subsmaples, from 2Mbvp to 12Mbp, of *Ameles dumonti* didn't output sufficient good assembly. Every BLAST of contigs from those analysis didn't match Mantodea mitochondrial genome. For this reason I follow this pipeline:
+ Assemble __Amedum__ with all trimmed reads avaible
+ Blast with [blastn](https://anaconda.org/channels/bioconda/packages/blast/files?version=2.14.1&type=) all database contigs assembly against *Ameles_andrea* COX1 
+ Chose best balsted contig from different analysis (with different reads subsample).
+ Mapping short reads with [Minimap2](https://github.com/lh3/minimap2/releases/tag/v2.26) on raw assembly contig and polishing with [Hypo](https://github.com/kensung-lab/hypo)   

Blast:
```bash
#BLASTN 2.14.1+
#Reference: Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb Miller (2000), "A greedy algorithm for aligning DNA sequences", J Comput Biol 2000; 7(1-2):203-14.
makeblastdb -in contigs.fasta -dbtype nucl -parse_seqids
blastn -query cox1ameand.fasta -db contigs.fasta
```

Mapping:
```bash
#script of Amedum with Minimap2
minimap2 -ax sr --MD -t 6 node_18_runall.fasta Amedum_1_paired.fastq  Amedum_2_paired.fastq > Amedum_raw_sr.sam
samtools view -Sb Amedum_raw_sr.sam > Amedum_raw_sr.bam
rm Amedum_raw_sr.sam
samtools sort -@6 -o Amedum_raw_sr_sorted.bam Amedum_raw_sr.bam
samtools index Amedum_raw_sr_sorted.bam
rm Amedum_raw_sr.bam

#run script
bash mapping_Amedum.sh
samtools flagstat Amedum_raw_sr_sorted.bam
```

Polishing: To run __Hypo__ it's necessary to obtain short reads coverage on assembly, using [mosdepth](https://anaconda.org/channels/bioconda/packages/mosdepth/overview)
```bash

#Installing and using mosdepth
conda install bioconda::mosdepth
mosdepth -n --fast-mode --by 500 Amedum_mito_coverage Amedum_raw_sr_sorted.bam
cat Amedum_mito_coverage.mosdepth.summary.txt 

#Installing and using Hypo
conda create -n hypo_env -c bioconda -c conda-forge hypo python=3.8
conda activate hypo_env
hypo -d node_18_runall.fasta -r Amedum_1_paired.fastq Amedum_2_paired.fastq -s 15678 -c 45 -b Amedum_raw_sr_sorted.bam -t 10 -o Amedum_mito_polished.fasta

#Annotate assembly with MITOS2 and then modify headers
sed 's/^>.*; *\([^;]*\) */>59T3_Amedum_[gene=\1]/' 59T3_Amedum_mitos.fasta > 59T3_Amedum.fasta
```

# Assembly results

| Method | Species | Code | Length(bp) | Coverage | Reads(bp) | 
| :---: | :---:	| :---: | :---: | :---: | :---:   
| lcWGS | *Ameles_assoi* | Ameass | 15634 | >10x | 1M |
| lcWGS | *Ameles_decolor* | Amedec | 15664 | >10x | 2M |
| lcWGS | *Apteromantis_aptera* | Aptapt | 16228 | >10x | 2M |
| lcWGS | *Pseudoyersinia_betancuriae* | Psebet | 15618 | >10x | 2M |
| SPAdes | *Ameles_spallanzania* | Amespa2 | 15714 | 12x | 2M |
| SPAdes | *Ameles_picteti* | Amepic | 15652 | 16x | 4M |
| SPAdes | *Rivetina_baetica* | Rivbae | 15613 | 16x | 4M |
| SPAdes | *Ameles_dumonti* | Amedum | 15678 | 8x | all |
