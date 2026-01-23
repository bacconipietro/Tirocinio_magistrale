# lcWGS MT Sweep Pipeline (Insect Mitogenomes)

This pipeline starts from raw short reads and performs:
1) trimming (Trimmomatic)
2) QC (FastQC)
3) MT assembly across a sweep of downsampled read counts (SPAdes)
4) MitoZ validation of mitochondrial contigs
5) Bowtie2 → samtools depth mean coverage calculation
6) automatic selection of the *minimum* sweep level achieving ≥10× mean depth (default)
7) TSV summary across all sweep runs
8) optional: remove mt reads, nuclear assembly, BUSCO (default insecta_odb10)

UCE/Phyluce is intentionally not enabled in the main environment to avoid conflicts with your existing Phyluce environment (`mosq-uce`).

---

## Installation

### Create the MT-sweep environment
```bash
mamba env create -f lcWGS_mt_sweep_env.yaml
mamba activate lcWGS_mt_sweep
```
