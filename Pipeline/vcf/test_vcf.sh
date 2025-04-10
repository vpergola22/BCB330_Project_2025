#!/bin/bash

python /ddn_exa/campbell/vpergola/Data/Pipeline/vcf/process_vcf2.py \
  --dir "/ddn_exa/campbell/vpergola/Data/Tumour/wes/song2022/processed2/variant_calling/deepvariant" \
  --output "/ddn_exa/campbell/vpergola/Data/Pipeline/vcf/song2022_vcf2.csv" \
  --chrom "CHROM" \
  --pos "POS" \
  --ref "REF" \
  --alt "ALT_1" \
  --samplesheet "/ddn_exa/campbell/vpergola/Data/Tests/id_handling/filtered_ids.csv"
