# LPS_typing_Illumina
Bioinformatics pipeline for *Pasteurella multocida* LPS typing using Illumina sequencing data

- [Quick start](#quick-start)
- [Overall pipeline](#overall-pipeline)
- [Step-by-step user guide](#step-by-step-user-guide)
- [Database setup](#database-setup)
- [Optional parameters](#optional-parameters)
- [Output files](#structure-of-the-output-folders)
- [Advanced use](#advanced-use)
- [Acknowledgements / citations / credits](#acknowledgements--citations--credits)

---

## Quick start

**Requirements:** [Nextflow](https://www.nextflow.io/docs/latest/install.html) ≥ 23.04, and either [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) or [Apptainer](https://apptainer.org/docs/user/latest/).

```bash
# 1. Clone the repository
git clone https://github.com/vmurigneu/LPS_typing_Illumina.git
cd LPS_typing_Illumina

# 2. Obtain large databases (see Database setup below)
#    The LPS and Kaptive databases are already included in databases/

# 3. Create your samplesheet (see Samplesheet section below)

# 4. Run (local execution with Apptainer)
nextflow run main.nf -profile apptainer \
  --samplesheet samplesheet/samples.csv \
  --outdir results

# On a SLURM cluster, add the slurm profile and your account name:
nextflow run main.nf -profile apptainer,slurm \
  --samplesheet samplesheet/samples.csv \
  --outdir results \
  --slurm_account YOUR_ACCOUNT
```

---

## Overall pipeline

### 1. Read trimming

Raw Illumina reads are trimmed using [fastp](https://github.com/OpenGene/fastp) v0.24.0.

### 2. Illumina reads quality metrics

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) v0.12.1 computes quality metrics for each sample on the trimmed reads. [MultiQC](https://github.com/MultiQC/MultiQC) v1.28 produces a combined report across all samples.

### 3. Genome assembly using Shovill

Paired-end reads are assembled using [Shovill](https://github.com/tseemann/shovill) v1.4.2, which uses SPAdes at its core.

### 4. Assembly quality assessment with QUAST

[QUAST](https://quast.sourceforge.net/quast.html) v5.2.0 computes assembly metrics on the polished assemblies.

### 5. Assembly quality assessment with CheckM

[CheckM](https://github.com/Ecogenomics/CheckM) v1.2.2 (command `lineage_wf`) estimates genome completeness and contamination from marker genes.

### 6. Kraken2/Bracken taxonomy classification

Trimmed reads are classified with [Kraken2](https://github.com/DerrickWood/kraken2) v2.1.3, and abundance is re-estimated at species level with [Bracken](https://github.com/jenniferlu717/Bracken) v3.0. The default database is the PlusPF index (Standard + RefSeq protozoa and fungi), see the [database index page](https://benlangmead.github.io/aws-indexes/k2) for details.

### 7. LPS typing using Kaptive

The LPS type is determined from the assembled genome using [Kaptive](https://kaptive.readthedocs.io/en/latest/) v3 with the 9-LPS database (included in `databases/kaptive3_LPS_db_v1/`).

### 8. Variant calling using Snippy

- Trimmed reads are mapped to the LPS reference sequence identified by Kaptive using BWA-mem as part of [Snippy](https://github.com/tseemann/snippy) v4.6.0.
- Snippy calls substitutions and indels relative to the reference LPS locus.
- [SnpEff](https://pcingola.github.io/SnpEff/) annotates variant effects on genes and proteins.
- High-impact variants (frameshift, stop_gained) are extracted to a separate file.

### 9. MLST typing

[mlst](https://github.com/tseemann/mlst) scans assemblies against the PubMLST scheme `pmultocida_2` (RIRDC) by default. The scheme can be changed with `--mlst_scheme`.

### 10. petG detection

[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) v2.17.0 searches the assembly for *petG* using the reference sequence `petG_X73_NZ_CM001580.fasta` (included in `databases/LPS/`). petG is reported as present when a hit spans > 1570 bp at ≥ 95% identity.

### 11. LPS subtype report

The pipeline generates a subtype report (`10_Illumina_subtype_report.tsv`) that matches Snippy variants against the [LPS subtype database](databases/LPS/LPS_subtype_database_v2.txt). A variant is reported when it matches the database at the same position, reference allele, and alternate allele in the same reference sequence.

When phenotype columns are present in the subtype database, the report assigns the corresponding LPS phenotype and its description from `phenotype_lookup.tsv`. The report also includes MLST sequence type and petG presence when those steps are run.

### 12. Genome annotation using Bakta

[Bakta](https://github.com/oschwengers/bakta) v1.12.0 annotates the genome assemblies. The default database is v6.0 (2025-02-24), available from [Zenodo record 10.5281/zenodo.14916843](https://zenodo.org/records/14916843).

### 13. Antimicrobial resistance genes

[AMRFinderPlus](https://github.com/ncbi/amr) v4.0.23 identifies AMR genes in the assemblies. The tested database version is 2025-03-25.1.

---

## Step-by-step user guide

### 1. Clone the repository

```bash
git clone https://github.com/vmurigneu/LPS_typing_Illumina.git
cd LPS_typing_Illumina
```

The repository includes three key files:

**a) `nextflow.config`** — default parameters and container profiles. Nextflow reads this automatically from the working directory. Select a container engine with `-profile singularity` or `-profile apptainer`. Add `-profile slurm` on SLURM clusters. Container images are pulled automatically and cached in `./singularity/` by default.

**b) `main.nf`** — the pipeline code. Not user-modifiable.

**c) `nextflow.sh`** — an optional SLURM submission script template. Edit it to set your account, partition, samplesheet path, and output directory, then submit with `sbatch nextflow.sh`.

### 2. Obtain databases

See the [Database setup](#database-setup) section below.

### 3. Prepare a samplesheet

The samplesheet is a comma-separated file listing each sample and its paired FASTQ files. Paths are resolved relative to the directory where you run `nextflow run` (the launch directory).

```bash
mkdir samplesheet
# Create/edit samplesheet/samples.csv:
```

```
sample_id,short_fastq_1,short_fastq_2
PM3034,fastq/PM3034_R1.fastq.gz,fastq/PM3034_R2.fastq.gz
PM3065,fastq/PM3065_R1.fastq.gz,fastq/PM3065_R2.fastq.gz
```

### 4. Run the pipeline

```bash
# Local execution with Apptainer:
nextflow run main.nf -profile apptainer \
  --samplesheet samplesheet/samples.csv \
  --outdir results

# SLURM cluster with Apptainer:
nextflow run main.nf -profile apptainer,slurm \
  --samplesheet samplesheet/samples.csv \
  --outdir results \
  --slurm_account YOUR_ACCOUNT

# Use -resume to restart from the last successful step:
nextflow run main.nf -profile apptainer,slurm \
  --samplesheet samplesheet/samples.csv \
  --outdir results \
  --slurm_account YOUR_ACCOUNT \
  -resume
```

---

## Database setup

The LPS reference data and Kaptive database are **included in the repository** under `databases/` and require no additional setup. Large third-party databases must be obtained separately; you can either download them automatically using pipeline flags, or provide pre-downloaded copies.

| Database | Required for | Default path | Tested version | Source | Approx. size | How to obtain |
|----------|-------------|--------------|---------------|--------|-------------|---------------|
| LPS references | Steps 7, 8, 10, 11 | `databases/LPS/` | Included in repo | This repository | < 1 MB | Included — no action needed |
| Kaptive3 LPS DB | Step 7 | `databases/kaptive3_LPS_db_v1/9lps.gbk` | Included in repo | This repository | < 200 KB | Included — no action needed |
| Kraken2 + Bracken | Step 6 | `databases/k2_pluspf_20240605` | k2_pluspf_20240605 | [AWS index page](https://benlangmead.github.io/aws-indexes/k2) | ~ 70 GB | Add `--download_kraken_db`, or download manually (see below) |
| CheckM | Step 5 | `databases/checkm_data_2015_01_16` | 2015_01_16 | [CheckM installation docs](https://github.com/Ecogenomics/CheckM/wiki/Installation) | ~ 1.4 GB | Add `--download_checkm_db`, or download manually (see below) |
| Bakta | Step 12 | `databases/bakta_db/db` | v6.0 (2025-02-24) | [Zenodo 10.5281/zenodo.14916843](https://zenodo.org/records/14916843) | ~ 65 GB | Add `--download_bakta_db`, or download manually (see below) |
| AMRFinderPlus | Step 13 | `databases/amrfinderplus/amrfinderplus_db/latest` | 2025-03-25.1 | [NCBI AMRFinderPlus wiki](https://github.com/ncbi/amr/wiki/AMRFinderPlus-database) | ~ 300 MB | Add `--download_amrfinder_db`, or download manually (see below) |

> **Reproducibility note:** The `--download_*` flags retrieve the latest available database version, which may differ from the tested versions listed above and could affect results. For reproducible analyses, use pinned database versions and specify their paths explicitly with `--kraken_db`, `--checkm_db`, `--bakta_db`, and `--amrfinder_db`.

### Downloading databases via pipeline flags

Add the relevant flag(s) on the first run. The database will be downloaded into `databases/` and reused automatically on subsequent runs (omit the flag after the first download).

```bash
# Download Kraken2/Bracken database (~70 GB, slow):
nextflow run main.nf -profile apptainer --samplesheet samplesheet/samples.csv \
  --outdir results --download_kraken_db

# Download CheckM database (~1.4 GB):
nextflow run main.nf -profile apptainer --samplesheet samplesheet/samples.csv \
  --outdir results --download_checkm_db

# Download Bakta database (~65 GB, slow):
nextflow run main.nf -profile apptainer --samplesheet samplesheet/samples.csv \
  --outdir results --download_bakta_db

# Download AMRFinderPlus database (~300 MB):
nextflow run main.nf -profile apptainer --samplesheet samplesheet/samples.csv \
  --outdir results --download_amrfinder_db
```

### Manual database download

If you prefer to download databases yourself (e.g. for speed or reproducibility), download them to any location and point the pipeline to them:

```bash
# Kraken2 (tested version k2_pluspf_20240605):
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20240605.tar.gz
tar -xvzf k2_pluspf_20240605.tar.gz -C /path/to/databases/k2_pluspf_20240605/
# Then run with: --kraken_db /path/to/databases/k2_pluspf_20240605

# CheckM (tested version checkm_data_2015_01_16):
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
mkdir -p /path/to/databases/checkm_data_2015_01_16
tar -xvzf checkm_data_2015_01_16.tar.gz -C /path/to/databases/checkm_data_2015_01_16/
# Then run with: --checkm_db /path/to/databases/checkm_data_2015_01_16

# Bakta (tested version v6.0 from Zenodo):
wget https://zenodo.org/records/14916843/files/db.tar.gz
tar -xvzf db.tar.gz -C /path/to/databases/bakta_db/
# Then run with: --bakta_db /path/to/databases/bakta_db/db

# AMRFinderPlus (tested version 2025-03-25.1):
# Use amrfinder_update within the container, or specify an existing database directory
# Then run with: --amrfinder_db /path/to/databases/amrfinderplus/2025-03-25.1
```

---

## Optional parameters

### 1. Read trimming

* `--skip_fastp`: skip the read trimming step (default=false). Not recommended.

### 2. FastQC reads quality metrics

* `--skip_fastqc`: skip the FastQC step (default=false)
* `--skip_summary_fastqc`: skip the MultiQC summary step (default=false)

### 3. Genome assembly

* `--shovill_threads`: number of threads for Shovill (default=4)
* `--shovill_args`: additional Shovill parameters (default="", example: `"--minlen 200 --mincov 10"`)
* `--genome_size`: estimated genome size (default="2.3M")

### 4. Assembly quality assessment with QUAST

* `--skip_quast`: skip the QUAST step (default=false)
* `--quast_threads`: number of threads for QUAST (default=2)

### 5. Assembly quality assessment with CheckM

* `--skip_checkm`: skip the CheckM step (default=false)
* `--checkm_db`: path to the CheckM database folder (default=`databases/checkm_data_2015_01_16`)
* `--download_checkm_db`: download the CheckM database automatically (default=false)

### 6. Kraken2/Bracken taxonomy classification

* `--skip_kraken`: skip the Kraken2/Bracken classification step (default=false)
* `--kraken_db`: path to the Kraken2 database folder (default=`databases/k2_pluspf_20240605`)
* `--download_kraken_db`: download the Kraken2 database automatically (default=false)

### 7. LPS typing using Kaptive

* `--skip_kaptive3`: skip the Kaptive typing step (default=false). Note: skipping this also skips variant calling (Snippy) and petG detection.
* `--kaptive_db_9lps`: path to the Kaptive database file (default=`databases/kaptive3_LPS_db_v1/9lps.gbk`)

### 8. Variant calling using Snippy

* `--skip_snippy`: skip the Snippy variant calling step (default=false)
* `--snippy_threads`: number of threads for Snippy (default=6)
* `--snippy_args`: additional Snippy parameters (default="")
* `--reference_LPS_directory`: path to the LPS reference directory containing `reference_LPS.txt`, `LPS_subtype_database_v2.txt`, `petG_X73_NZ_CM001580.fasta`, and FASTA/GB files for all LPS types (default=`databases/LPS`)

### 9. MLST typing

* `--skip_mlst`: skip the MLST typing step (default=false)
* `--mlst_scheme`: MLST typing scheme (default="pmultocida_2")

### 10. petG detection

* `--skip_petg`: skip petG detection with BLAST (default=false)
* `--petg_threads`: number of threads for the BLAST step (default=2)
* `--petg_min_length`: minimum genomic hit span in bp for petG presence; hits must exceed this value (default=1570)
* `--petg_min_identity`: minimum percent identity for petG presence (default=95)

### 11. Report

The subtype report uses `LPS_subtype_database_v2.txt` and, when present, `phenotype_lookup.tsv` from `--reference_LPS_directory`. No additional parameters needed.

### 12. Genome annotation using Bakta

* `--skip_bakta`: skip genome annotation (default=false)
* `--bakta_threads`: number of threads for Bakta (default=8)
* `--bakta_db`: path to the Bakta database folder (default=`databases/bakta_db/db`)
* `--download_bakta_db`: download the latest Bakta database automatically (default=false)
* `--bakta_args`: additional Bakta parameters (default includes `--proteins databases/LPS/NC_002663_LPS.gb` for *Pasteurella multocida*)

### 13. AMR gene identification using AMRFinderPlus

* `--skip_amrfinder`: skip AMR gene identification (default=false)
* `--amrfinder_db`: path to the AMRFinderPlus database folder (default=`databases/amrfinderplus/amrfinderplus_db/latest`)
* `--download_amrfinder_db`: download the latest AMRFinderPlus database automatically (default=false)
* `--amrfinder_args`: additional AMRFinderPlus parameters (default="")

---

## Structure of the output folders

The pipeline creates a folder per sample (named by `sample_id`) inside `--outdir`, plus a combined `10_report` folder.

Each sample folder contains:

* **1_trimming:** Trimmed paired-end FASTQ files (`sample_id_R1_trimmed.fastq.gz`, `sample_id_R2_trimmed.fastq.gz`).
* **2_fastqc:** FastQC results — HTML report and ZIP archive for each read set.
* **3_assembly:** Shovill assembly output. See [Shovill output docs](https://github.com/tseemann/shovill?tab=readme-ov-file#output-files).
  * Final assembly in FASTA format (`sample_id_contigs.fa`)
  * SPAdes assembly graph (`sample_id_contigs.gfa`)
* **4_quast:** QUAST report (`sample_id_report.tsv`).
* **5_checkm:** CheckM lineage workflow results (`sample_id_checkm_lineage_wf_results.tsv`).
* **6_kraken:** Kraken2/Bracken taxonomy classification results. See output format details [here (Kraken2)](https://github.com/DerrickWood/kraken2/wiki/Manual#output-formats) and [here (Bracken)](https://ccb.jhu.edu/software/bracken/index.shtml?t=manual#format).
  * Kraken2 classification report (`sample_id_kraken2_report.txt`)
  * Kraken2 per-read assignments (`sample_id_kraken2.tsv.gz`)
  * Bracken species abundance (`sample_id_bracken_species.tsv`)
  * Bracken Kraken-style report (`sample_id_bracken_report.txt`)
* **7_kaptive_v3:** Kaptive output files. See [Kaptive output docs](https://kaptive.readthedocs.io/en/latest/Outputs.html).
  * LPS type results (`sample_id_kaptive_results.tsv`)
  * LPS sequence in FASTA format (`sample_id_kaptive_results.fna`)
* **8_snippy:** Snippy mapping and variant calling results. See [Snippy output docs](https://github.com/tseemann/snippy?tab=readme-ov-file#output-files).
  * BAM alignment file (`sample_id_snps.bam_mapped.bam` and `.bai` index)
  * Unfiltered variants VCF (`sample_id_clair_snps.raw.vcf`)
  * Filtered variants VCF (`sample_id_snps.filt.vcf`)
  * Variant summary in tabular format (`sample_id_snps.tab`)
  * High-impact variants (frameshift, stop_gained) (`sample_id_snps.high_impact.tab`)
* **9_mlst:** MLST typing result (`sample_id_mlst_pmultocida_rirdc.csv`).
* **13_petG:** petG BLAST output files.
  * BLAST tabular output for all hits (`sample_id_petG_blast.tsv`)
  * BLAST tabular output for accepted hits (`sample_id_petG_blast.filtered.tsv`)
  * FASTA sequences for accepted hits (`sample_id_petG_hits.fasta`)
  * petG presence summary (`sample_id_petG_summary.tsv`)
* **11_bakta:** Bakta annotation output. See [Bakta output docs](https://github.com/oschwengers/bakta?tab=readme-ov-file#output).
  * Annotations in GenBank format (`sample_id_bakta.gbff`)
  * Inference metrics (`sample_id_bakta.inference.tsv`)
  * Annotation summary (`sample_id_bakta.txt`)
* **12_amrfinder:** AMRFinderPlus results (`sample_id_amrfinder.tsv`). See [output format](https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus#output-format).

The `10_report` folder contains combined results across all samples:

* MultiQC report (`2_Illumina_multiqc_report.html`) and general statistics (`2_Illumina_multiqc_general_stats.txt`)
* Shovill assembly statistics: coverage, contig count, assembly size (`3_Illumina_shovill_stats.tsv`)
* Combined QUAST report (`4_Illumina_quast_report.tsv`)
* Combined CheckM results (`5_Illumina_checkm_lineage_wf_results.tsv`)
* Kraken/Bracken taxonomy results:
  * *P. multocida* abundance (`6_Illumina_bracken_pasteurella_multocida_species_abundance.tsv`)
  * Most abundant species (`6_Illumina_bracken_most_abundant_species.tsv`)
* Kaptive results (`7_Illumina_kaptive_results.tsv`)
* Snippy variant results:
  * All variants (`8_Illumina_snippy_snps.tsv`)
  * High-impact variants only (`8_Illumina_snippy_snps.high_impact.tsv`)
* MLST results (`9_Illumina_mlst.csv`)
* Subtype report (`10_Illumina_subtype_report.tsv`). Column descriptions:
  * **SAMPLE**: sample identifier
  * **MLST**: MLST sequence type
  * **TYPE**: LPS type assigned by Kaptive (or `untypeable`)
  * **SUBTYPE**: LPS subtype from the subtype database
  * **VARTYPE**: description of the variant
  * **ISOLATE_DATABASE**: reference isolate from the subtype database
  * **CHROM**: reference sequence name for the LPS locus
  * **POS**: variant position in the reference
  * **REF**: reference allele
  * **ALT**: alternate allele identified in the sample
  * **GENE**: gene containing the variant
  * **PHENOTYPE**: LPS phenotype from the subtype database, when available
  * **PHENOTYPE_DESCRIPTION**: phenotype description from `phenotype_lookup.tsv`, when available
  * **PETG_PRESENT**: petG presence (`yes` when present, blank otherwise)
  * **NOTE**: subtype database note for the matched variant, when available
* AMRFinderPlus combined results (`12_Illumina_amrfinder.tsv`)

---

## Advanced use

### Running in assembly-only mode for other organisms

The default parameters are optimised for *Pasteurella multocida*. LPS typing, variant calling, and petG detection are species-specific. To use the pipeline for genome assembly, QC, taxonomy, and MLST for another species, skip the *P. multocida*-specific steps:

```bash
nextflow run main.nf -profile apptainer \
  --samplesheet samplesheet/samples.csv \
  --outdir results \
  --genome_size 2.5M \
  --mlst_scheme your_scheme \
  --skip_kaptive3 \
  --skip_bakta \
  --skip_amrfinder
```

Note: `--skip_kaptive3` automatically skips Snippy variant calling and petG detection.

### Using a custom cluster configuration

For clusters with non-standard SLURM partitions, memory limits, or other requirements, create a custom config file and pass it with `-c`:

```bash
# my_cluster.config
process {
    executor = 'slurm'
    clusterOptions = '--account=my_account --partition=high_mem'
    time = '12h'
    withLabel: high_memory { memory = 512.GB }
}
```

```bash
nextflow run main.nf -profile apptainer \
  --samplesheet samplesheet/samples.csv \
  --outdir results \
  -c my_cluster.config
```

---

## Acknowledgements / citations / credits

Please cite the following tools when using this pipeline:

- [fastp](https://github.com/OpenGene/fastp)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](https://github.com/MultiQC/MultiQC)
- [Shovill](https://github.com/tseemann/shovill) / [SPAdes](https://github.com/ablab/spades)
- [QUAST](https://quast.sourceforge.net/quast.html)
- [CheckM](https://github.com/Ecogenomics/CheckM)
- [Kraken2](https://github.com/DerrickWood/kraken2)
- [Bracken](https://github.com/jenniferlu717/Bracken)
- [Kaptive](https://kaptive.readthedocs.io/en/latest/)
- [Snippy](https://github.com/tseemann/snippy) / [SnpEff](https://pcingola.github.io/SnpEff/)
- [mlst](https://github.com/tseemann/mlst)
- [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- [Bakta](https://github.com/oschwengers/bakta)
- [AMRFinderPlus](https://github.com/ncbi/amr)

Pipeline developed by Valentine Murigneux and Julian Zaugg.
