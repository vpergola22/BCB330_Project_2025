#!/bin/bash

nextflow run nf-core/oncoanalyser  \
    -profile singularity  \
    -revision 1.0.0  \
    --mode wgts  \
    --genome GRCh38_hmf \
    --input /ddn_exa/campbell/vpergola/github/nfcore_runs/WES/song2022/oncoanalyzer_samplesheet.csv \
    --outdir /ddn_exa/campbell/vpergola/Data/Tumour/wes/song2022/processed