#!/bin/bash
#SBATCH --partition=long
#SBATCH --output=job_run_rnaseq_denovo_Ecry_%j.out
#SBATCH --error=job_run_rnaseq_denovo_Ecry_%j.err
#SBATCH --mem-per-cpu=5G
#SBATCH --cpus-per-task=1

source activate rnaseq_denovo_env

snakemake --profile workflow/profiles/ --rerun-incomplete 

echo Complete!

