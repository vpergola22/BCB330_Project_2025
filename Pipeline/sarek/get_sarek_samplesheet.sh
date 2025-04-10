#!/bin/bash

python /ddn_exa/campbell/vpergola/Data/Tests/get_sarek_samplesheet.py \
  --input "/ddn_exa/campbell/vpergola/Data/Tumour/wes/song2022/metadata/PRJNA754592.runinfo_ftp.tsv" \
  --output "/ddn_exa/campbell/vpergola/Data/Pipeline/song2022_sarek_samplesheet.csv" \
  --filter "/ddn_exa/campbell/vpergola/Data/Tests/id_handling/filtered_ids.csv"
