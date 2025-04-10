#!/bin/bash

nextflow run nf-core/fetchngs \
   -profile conda \
   --input "/ddn_exa/campbell/vpergola/github/nfcore_runs/WES/song2022/ids.csv" \
   --outdir "/ddn_exa/campbell/vpergola/Data/Tumour/wes/song2022_ver2" \
   --nf_core_pipeline rnaseq



