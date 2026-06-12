#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --job-name=LPS_pipeline
#SBATCH --output=./%j_LPS_pipeline.out
#SBATCH --error=./%j_LPS_pipeline.err
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1

# Edit the paths below before submitting
SAMPLESHEET=/path/to/samplesheet/samples.csv
OUTDIR=/path/to/results

# Run the pipeline (-resume restarts from the last successful step if the job is interrupted)
nextflow run main.nf \
  -profile apptainer,slurm \
  --samplesheet ${SAMPLESHEET} \
  --outdir ${OUTDIR} \
  --slurm_account YOUR_ACCOUNT \
  -resume
