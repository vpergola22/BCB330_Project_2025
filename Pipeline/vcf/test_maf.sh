#!/bin/bash

python /ddn_exa/campbell/vpergola/Data/Pipeline/vcf/process_maf.py \
  --dir "/ddn_exa/campbell/vpergola/Data/Pipeline/meskit/final_mafs" \
  --output "/ddn_exa/campbell/vpergola/Data/Joined/song2022_maf.csv" \
  --samplesheet "/ddn_exa/campbell/vpergola/Data/Pipeline/ids/id_handling/filtered_ids.csv"
