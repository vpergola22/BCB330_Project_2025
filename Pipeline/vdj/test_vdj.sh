#!/bin/bash

python /ddn_exa/campbell/vpergola/Data/Tcell/sc/process_vdj.py /ddn_exa/campbell/vpergola/Data/Tcell/sc/song2022 /ddn_exa/campbell/vpergola/Data/Tests/song2022_vdj.csv \
  --chain_col "chain" \
  --v_col "v_gene" \
  --j_col "j_gene" \
  --barcode_col "barcode" \
  --contig_col "contig_id"
