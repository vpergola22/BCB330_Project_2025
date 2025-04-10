#!/bin/bash

# Set input and output directories
VCF_INPUT_DIR="/ddn_exa/campbell/vpergola/Data/Pipeline/meskit/vcf_inputs/run_1"
MAF_OUTPUT_DIR="/ddn_exa/campbell/vpergola/Data/Pipeline/meskit/meskit_input"

# Reference and VEP paths
REF_FASTA="/ddn_exa/campbell/vpergola/.vep/homo_sapiens/112_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
VEP_PATH="/ddn_exa/campbell/vpergola/miniforge3/envs/vep/bin"

# Tumor ID list - IMPORTANT: Order must match the order of VCF files when sorted
TUMOR_IDS=(
    #SAMN20792896_SRX11758678 
    SAMN20792920_SRX11758680
    #SAMN20792851_SRX11758698
    #SAMN20792930_SRX11758708
    #SAMN20792935_SRX11758710
    #SAMN20792938_SRX11758711
    #SAMN20792942_SRX11758712
    #SAMN20792945_SRX11758713
    # Add more tumor IDs as needed
    # IMPORTANT make sure no quotes or commas are used, had to do them one by one 
    # TODO find out how to separate items in a list in bash
)

# Create output directory if it doesn't exist
mkdir -p "$MAF_OUTPUT_DIR"

# Log file to track processing
LOG_FILE="${MAF_OUTPUT_DIR}/vcf2maf_conversion.log"

# Clear previous log
> "$LOG_FILE"

# Process each VCF file in the input directory
# Use a counter to match tumor IDs
i=0
for vcf_file in "$VCF_INPUT_DIR"/SRX*.deepvariant.vcf; do
    # Extract filename without path
    filename=$(basename "$vcf_file")
    
    # Create output MAF filename
    maf_file="${MAF_OUTPUT_DIR}/${filename%.deepvariant.vcf}.vep.maf"

    # Check if we have a tumor ID for this file
    if [ $i -lt ${#TUMOR_IDS[@]} ]; then
        tumor_id="${TUMOR_IDS[$i]}"
    else
        echo "ERROR: No tumor ID available for $filename" | tee -a "$LOG_FILE"
        continue
    fi
    
    # Run vcf2maf conversion
    echo "Processing $filename ..." | tee -a "$LOG_FILE"
    
    perl /ddn_exa/campbell/vpergola/.vep/homo_sapiens/112_GRCh38/mskcc-vcf2maf-f6d0c40/vcf2maf.pl \
        --input-vcf "$vcf_file" \
        --output-maf "$maf_file" \
        --ref-fasta "$REF_FASTA" \
        --vep-path "$VEP_PATH" \
        --ncbi-build GRCh38 \
        --tumor-id "$tumor_id" \
        --retain-info "AD,VAF" \
        --retain-fmt "AD,VAF" \
        --retain-ann "AD_REF,AD_ALT" \
        --verbose \
        2>&1 | tee -a "$LOG_FILE"
    
    # Check conversion status
    if [ $? -eq 0 ]; then
        echo "Successfully converted $filename to MAF" | tee -a "$LOG_FILE"
    else
        echo "ERROR: Failed to convert $filename" | tee -a "$LOG_FILE"
    fi

    # Increment counter
    ((i++))

done

# Summary
echo "Batch VCF to MAF conversion completed. Check log at $LOG_FILE" | tee -a "$LOG_FILE"