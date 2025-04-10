#!/bin/bash
#SBATCH --job-name=nf-sarek-scrnaseq
#SBATCH --output=nextflow_pipeline_%j.log
#SBATCH --error=nextflow_pipeline_%j.err
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=256GB

# Debugging settings
export APPTAINER_TMPDIR=/ddn_exa/campbell/vpergola/github/nfcore_extras/scratch
export APPTAINER_CACHEDIR=/ddn_exa/campbell/vpergola/github/nfcore_extras/cache

# config
config="/ddn_exa/campbell/vpergola/nf_configs/galen.config"

nextflow run nf-core/sarek \
    -r 3.5.0 \
    -profile singularity \
    --input /ddn_exa/campbell/vpergola/Data/Pipeline/song2022_sarek_samplesheet.csv \
    --outdir /ddn_exa/campbell/vpergola/Data/Tumour/wes/song2022/processed2 \
    --genome GATK.GRCh38 \
    --tools deepvariant \
    -process.cache lenient \
    -resume \
    -c $config \
    --reference /ddn_exa/campbell/vpergola/github/nfcore_runs/WES/song2022/work/f9/1227745557ca4d1f9620a05cb2e3bf/Homo_sapiens_assembly38.fasta  \
    -executor.queueSize 20