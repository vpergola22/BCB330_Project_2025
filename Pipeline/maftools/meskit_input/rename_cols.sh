#!/bin/bash

MAF_INPUT_DIR="/ddn_exa/campbell/vpergola/Data/Pipeline/meskit/meskit_input/not_done"
MAF_OUTPUT_DIR="/ddn_exa/campbell/vpergola/Data/Pipeline/meskit/final_mafs/"

mkdir -p "$MAF_OUTPUT_DIR"

for maf_file in "$MAF_INPUT_DIR"/*.vep.maf; do
    filename=$(basename "$maf_file")
    output_file="${MAF_OUTPUT_DIR}/${filename}"
    
    # Use awk to process the MAF file
    awk 'BEGIN {FS=OFS="\t"}
    NR==2 {
        # Rename first VAF column to VAF_empty
        for(i=1;i<=NF;i++) {
            if($i=="VAF") {
                $i="VAF_empty"
                break
            }
        }

        # Find column indices
        for(i=1;i<=NF;i++) {
            if($i=="t_VAF") vaf_col=i
            if($i=="t_AD") ad_col=i
        }
        
        # Rename columns
        $vaf_col="VAF"
        
        # Add new columns for ref and alt allele depths
        $0 = $0 "\tRef_allele_depth\tAlt_allele_depth"
        print
        next
    }
    
    {
        # Split AD column
        split($ad_col, ad_vals, ",")
        
        # Replace VAF column
        $vaf_col=$vaf_col
        
        # Add ref and alt allele depths
        $0 = $0 "\t" ad_vals[1] "\t" ad_vals[2]
        
        print
    }' "$maf_file" > "$output_file"
    
    echo "Processed $filename"
done

echo "MAF file processing complete."