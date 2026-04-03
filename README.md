# MBGAP — Mykaella's Bacterial Genome Analysis Pipeline

**Version 3.0**  
Developed by **Jean Phellipe Marques do Nascimento**  
Laboratório de Vigilância Genômica — LACEN-AL

---

## What is MBGAP?

MBGAP is an automated pipeline for bacterial genome analysis from Illumina **paired-end** next-generation sequencing (NGS) data. With a single command, it runs the following steps in an integrated manner:

1. **Raw read quality control** (FastQC)
2. **Adapter and quality trimming** (Trimmomatic) — optional
3. **Genome assembly** — SPAdes, Shovill, or Unicycler
4. **Assembly quality assessment** (QUAST + assembly-scan)
5. **Taxonomic identification** (GAMBIT and/or GTDB-Tk)
6. **Genome annotation** (Prokka + Bakta)
7. **Antimicrobial resistance gene detection** (AMRFinder)
8. **Plasmid identification** (PlasmidFinder)
9. **Consolidated quality report** (MultiQC)

---

## Before You Start — Prerequisites

### Operating System

MBGAP is designed to run on **Linux** (Ubuntu/Debian recommended). It is not directly compatible with Windows (you may use WSL2 on Windows if needed).

### Conda

All tools are managed via **Conda**. If you do not have Conda installed:

1. Download the Miniconda installer:
   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   ```
2. Restart your terminal or run `source ~/.bashrc`.

### Required Conda Environments

MBGAP uses **two Conda environments**:

| Environment | Included Tools |
|---|---|
| `bioinfo` (default) | FastQC, Trimmomatic, SPAdes, Shovill, Unicycler, QUAST, Prokka, Bakta, AMRFinder, PlasmidFinder, GAMBIT, assembly-scan, MultiQC |
| `gtdbtk` (default) | GTDB-Tk |

> GTDB-Tk requires a separate environment due to dependency conflicts with other tools.

You can use different names for the environments — just provide them via `--conda-env` and `--gtdbtk-env`.

**Suggested installation of the main environment** (adjust as needed):
```bash
conda create -n bioinfo -c conda-forge -c bioconda \
  fastqc trimmomatic spades shovill unicycler quast \
  prokka bakta amrfinder plasmidfinder gambit assembly-scan multiqc
```

### Required Databases

You must download and configure **three databases** before using MBGAP:

| Database | Tool | How to obtain |
|---|---|---|
| GAMBIT DB | GAMBIT | [github.com/jlumpe/gambit](https://github.com/jlumpe/gambit) |
| Bakta DB | Bakta | `bakta_db download --output /path/to/bakta_db` |
| PlasmidFinder DB | PlasmidFinder | [bitbucket.org/genomicepidemiology/plasmidfinder_db](https://bitbucket.org/genomicepidemiology/plasmidfinder_db) |

---

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/YOUR_USERNAME/MBGAP.git
   cd MBGAP
   ```
   Or download the `MBGAP_v3_0.sh` file directly.

2. Make the script executable:
   ```bash
   chmod +x MBGAP_v3_0.sh
   ```

3. That's it! No additional installation is needed beyond the Conda environments described above.

---

## How to Use

### Input Modes

You can provide data in two ways:

**Option A — Directory of FASTQ files**

If your files follow the naming pattern `SAMPLE_R1_001.fastq.gz` / `SAMPLE_R2_001.fastq.gz` (or `SAMPLE_1.fastq.gz` / `SAMPLE_2.fastq.gz`):

```bash
bash MBGAP_v3_0.sh \
  --input /path/to/fastqs/ \
  --output /path/to/results/ \
  --assembler spades \
  --gambit-db /path/to/gambit_db/ \
  --bakta-db /path/to/bakta_db/ \
  --plasmidfinder-db /path/to/plasmidfinder_db/
```

**Option B — CSV Samplesheet (recommended for multiple samples with metadata)**

Create a CSV file with information for each sample:

```csv
sample,r1,r2,genus,species,gram,gsize
ISO001,/data/ISO001_R1.fastq.gz,/data/ISO001_R2.fastq.gz,Escherichia,coli,-,5.2M
ISO002,/data/ISO002_R1.fastq.gz,/data/ISO002_R2.fastq.gz,Klebsiella,pneumoniae,-,5.5M
ISO003,/data/ISO003_R1.fastq.gz,/data/ISO003_R2.fastq.gz,,,,4.9M
```

> **Required fields:** `sample`, `r1`, `r2`  
> **Optional fields:** `genus`, `species`, `gram` (`+`, `-`, or `?`), `gsize` (required for Shovill)  
> Leave optional fields blank when unknown — GAMBIT can automatically infer taxonomy.

```bash
bash MBGAP_v3_0.sh \
  --samplesheet /path/to/samples.csv \
  --output /path/to/results/ \
  --assembler spades \
  --gambit-db /path/to/gambit_db/ \
  --bakta-db /path/to/bakta_db/ \
  --plasmidfinder-db /path/to/plasmidfinder_db/
```

---

## Full Parameter Reference

### Input and Output

| Parameter | Description |
|---|---|
| `--input DIR` | Directory containing paired-end FASTQ files |
| `--samplesheet FILE` | CSV file with per-sample metadata (replaces `--input`) |
| `--output DIR` | **Required.** Directory where results will be saved |

### Databases (all required)

| Parameter | Description |
|---|---|
| `--gambit-db DIR` | Path to the GAMBIT database |
| `--bakta-db DIR` | Path to the Bakta database |
| `--plasmidfinder-db DIR` | Path to the PlasmidFinder database |

### Assembly

| Parameter | Values | Default | Description |
|---|---|---|---|
| `--assembler` | `spades`, `shovill`, `unicycler` | `spades` | Assembler to use |
| `--gsize` | e.g. `5.0M`, `4500k` | — | Estimated genome size (required for Shovill when not provided in the CSV) |
| `--memory-spades` | integer (GB) | `50` | Maximum memory for SPAdes |
| `--ram-shovill` | integer (GB) | `50` | Maximum RAM for Shovill |

### Quality and Trimming

| Parameter | Values | Default | Description |
|---|---|---|---|
| `--trim` | `yes`, `no` | `no` | Enable Trimmomatic trimming |
| `--adapters FILE` | path to file | — | Trimmomatic adapter file (**required** if `--trim yes`) |

### Taxonomic Information (optional global defaults)

| Parameter | Values | Default | Description |
|---|---|---|---|
| `--genus` | text | — | Global bacterial genus for all samples |
| `--species` | text | — | Global bacterial species for all samples |
| `--gram` | `+`, `-`, `?` | `?` | Gram staining type |
| `--genetic-code` | integer | `11` | Genetic code for Prokka (11 = bacteria) |
| `--use-prokka-genus` | `yes`, `no` | `no` | Enable `--usegenus` flag in Prokka |
| `--bakta-complete` | `yes`, `no` | `no` | Inform Bakta that the genome is complete |

### Optional Tools

| Parameter | Values | Default | Description |
|---|---|---|---|
| `--run-gambit` | `yes`, `no` | `yes` | Run GAMBIT taxonomic identification |
| `--run-gtdbtk` | `yes`, `no` | `yes` | Run GTDB-Tk phylogenomic classification |
| `--run-multiqc` | `yes`, `no` | `yes` | Generate MultiQC report |

### Conda Environments

| Parameter | Default | Description |
|---|---|---|
| `--conda-env` | `bioinfo` | Name of the main Conda environment |
| `--gtdbtk-env` | `gtdbtk` | Name of the GTDB-Tk Conda environment |

### Threads

| Parameter | Default | Description |
|---|---|---|
| `--threads INT` | — | Set all threads at once |
| `--threads-fastqc INT` | `8` | Threads for FastQC |
| `--threads-trimmomatic INT` | `8` | Threads for Trimmomatic |
| `--threads-spades INT` | `10` | Threads for SPAdes |
| `--threads-shovill INT` | `10` | Threads for Shovill |
| `--threads-unicycler INT` | `10` | Threads for Unicycler |
| `--threads-quast INT` | `8` | Threads for QUAST |
| `--threads-prokka INT` | `8` | Threads for Prokka |
| `--threads-bakta INT` | `12` | Threads for Bakta |
| `--threads-amrfinder INT` | `8` | Threads for AMRFinder |
| `--threads-gtdbtk INT` | `10` | Threads for GTDB-Tk |

---

## Output Structure

After execution, the output directory (`--output`) will have the following structure:

```
results/
├── pipeline_log.txt                        # Full execution log
├── sample_taxonomy_resolution.tsv          # Per-sample taxonomy summary
├── taxonomic_identification_gambit.tsv     # GAMBIT results
├── taxonomic_identification_gtdbtk/        # GTDB-Tk results
├── multiqc_output/                         # Consolidated HTML report (MultiQC)
├── Assembly_final/                         # Final FASTA files for all samples
│   ├── ISO001_spades.fasta
│   └── ISO002_spades.fasta
├── ISO001/                                 # Per-sample directory
│   ├── ISO001_fastqc/                      # FastQC on raw reads
│   ├── ISO001_fastqc_trimmed/              # FastQC after trimming (if --trim yes)
│   ├── ISO001_trimmomatic/                 # Trimmed reads
│   ├── ISO001_spades/                      # Raw SPAdes output
│   ├── ISO001_quast/                       # Assembly statistics
│   ├── ISO001_spades.tsv                   # assembly-scan summary
│   ├── ISO001_prokka/                      # Prokka annotation
│   ├── ISO001_bakta/                       # Bakta annotation
│   ├── ISO001_amrfinder/                   # AMR genes (AMRFinder)
│   └── ISO001_plasmidfinder/               # Identified plasmids
└── ISO002/
    └── ...
```

---

## Usage Examples

### Example 1 — Basic run (no trimming, SPAdes, all defaults)

```bash
bash MBGAP_v3_0.sh \
  --input /data/fastqs/ \
  --output /results/run_01/ \
  --assembler spades \
  --gambit-db /databases/gambit/ \
  --bakta-db /databases/bakta/ \
  --plasmidfinder-db /databases/plasmidfinder/
```

### Example 2 — With trimming and global genus/species

```bash
bash MBGAP_v3_0.sh \
  --input /data/fastqs/ \
  --output /results/run_02/ \
  --assembler spades \
  --trim yes \
  --adapters /databases/TruSeq3-PE.fa \
  --genus Staphylococcus \
  --species aureus \
  --gram + \
  --gambit-db /databases/gambit/ \
  --bakta-db /databases/bakta/ \
  --plasmidfinder-db /databases/plasmidfinder/
```

### Example 3 — Samplesheet with Shovill and 20 threads

```bash
bash MBGAP_v3_0.sh \
  --samplesheet /data/my_samples.csv \
  --output /results/run_03/ \
  --assembler shovill \
  --threads 20 \
  --gambit-db /databases/gambit/ \
  --bakta-db /databases/bakta/ \
  --plasmidfinder-db /databases/plasmidfinder/
```

### Example 4 — Without GTDB-Tk (faster, no GTDB reference data needed)

```bash
bash MBGAP_v3_0.sh \
  --samplesheet /data/samples.csv \
  --output /results/run_04/ \
  --assembler spades \
  --run-gtdbtk no \
  --gambit-db /databases/gambit/ \
  --bakta-db /databases/bakta/ \
  --plasmidfinder-db /databases/plasmidfinder/
```

---

## Automatic Checkpoint and Resume

MBGAP saves the progress of each step. If execution is interrupted (power failure, error, etc.), simply run the same command again — already completed steps will be skipped automatically.

Checkpoint files are stored in `OUTPUT_DIR/.checkpoints/`.

---

## Troubleshooting

**Error: `Dependency not found in PATH: fastqc`**  
→ The Conda environment is not active or the tool was not installed. Check with `conda activate bioinfo` and then `which fastqc`.

**Error: `Could not activate conda environment: bioinfo`**  
→ Your environment has a different name. Check with `conda env list` and use `--conda-env YOUR_ENV_NAME`.

**Error: `gsize missing for sample X using Shovill`**  
→ When using `--assembler shovill`, the estimated genome size is required. Add the `gsize` column to your CSV or use `--gsize 5.0M`.

**Error: `--adapters is required when --trim yes`**  
→ Provide the path to your Illumina adapter file with `--adapters /path/to/TruSeq3-PE.fa`.

**Error: `Duplicate sample detected`**  
→ Two FASTQ files produced the same sample name, or there are duplicate rows in the CSV.

**The pipeline is running slowly**  
→ Increase the number of threads with `--threads 16` (or more, depending on your server) and, if needed, increase `--memory-spades` and `--ram-shovill`.

---

## Citation

If you use MBGAP in a scientific publication, please cite this repository and the individual tools used:

- **SPAdes:** Bankevich et al., 2012 — *J Comput Biol*
- **Shovill:** [github.com/tseemann/shovill](https://github.com/tseemann/shovill)
- **Unicycler:** Wick et al., 2017 — *PLOS Computational Biology*
- **FastQC:** Andrews, S. — [bioinformatics.babraham.ac.uk](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- **Trimmomatic:** Bolger et al., 2014 — *Bioinformatics*
- **QUAST:** Gurevich et al., 2013 — *Bioinformatics*
- **Prokka:** Seemann, 2014 — *Bioinformatics*
- **Bakta:** Schwengers et al., 2021 — *Microbial Genomics*
- **AMRFinder:** Feldgarden et al., 2021 — *Scientific Reports*
- **PlasmidFinder:** Carattoli et al., 2014 — *Antimicrobial Agents and Chemotherapy*
- **GAMBIT:** Lumpe & Bhatt — [github.com/jlumpe/gambit](https://github.com/jlumpe/gambit)
- **GTDB-Tk:** Chaumeil et al., 2019 — *Bioinformatics*
- **MultiQC:** Ewels et al., 2016 — *Bioinformatics*

---

## License

This project is licensed under the [MIT License](./LICENSE).  
Developed at the **Laboratório de Vigilância Genômica — LACEN-AL**.

---

*Questions or suggestions? Open an [issue](../../issues) in this repository.*
