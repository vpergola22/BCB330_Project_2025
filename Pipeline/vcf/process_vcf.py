import allel as allel
import pandas as pd
import os
import argparse
import re

'''
# original code
df = allel.vcf_to_dataframe("/ddn_exa/campbell/vpergola/Data/Reference/annotated_test.vcf") #TODO change the location and name of this

# FILTERING DATA TODO this may be different for different studies
# filter by filter_pass and check if mutation is blank
df = df[df['FILTER_PASS'] == True]

# df = df[df['REF'] != df['ALT_1']] not sure if this is needed CHECK 
df = df[df['ALT_1'] != "<*>"]

# get rid of columns we dont need
df = df.drop(columns=['ALT_2', 'ALT_3'])
'''

def process_vcf(directory, output_file, column_map):
    """
    Processes multiple VCF files in subdirectories and generates a binary presence-absence matrix.
    
    Parameters:
        directory (str): Directory where SRX subdirectories are located.
        output_file (str): Output CSV file path.
        column_map (dict): Dictionary containing column names for chrom, pos, ref, alt, and other necessary columns.
    """
    
    # Initialize an empty list to store VCF file paths
    vcf_files = []

    # Recursively search for VCF files in SRX subdirectories
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith("deepvariant.vcf.gz"): #TODO change after testing to use the vcf.gz file instead
                # Check if the file is within an SRX subdirectory
                if "deepvariant" in root:  # Ensure it's within the 'deepvariant' subdirectory
                    vcf_files.append(os.path.join(root, file))

    print(f"Found VCF files: {vcf_files}")  # Debugging: List the found VCF files

    master_df = pd.DataFrame()

    for vcf_file in vcf_files:
        # Extract sample ID robustly
        match = re.search(r"(SRX\d+)", vcf_file)
        if match:
            sample_id = match.group(1)
        else:
            print(f"Warning: Could not extract Sample_ID from {vcf_file}")
            continue  # Skip this file if sample ID is not found

        # Convert VCF to DataFrame
        try:
            df = allel.vcf_to_dataframe(vcf_file)
            print(df.head())  # Debugging
        except Exception as e:
            print(f"Error processing {vcf_file}: {e}")
            continue

        # Rename columns dynamically based on user input
        df = df.rename(columns={
            column_map['chrom_col']: "CHROM",
            column_map['pos_col']: "POS",
            column_map['ref_col']: "REF",
            column_map['alt_col']: "ALT"
        })

        # Apply filters (Modify as needed for other datasets)
        # TODO might have to modify this if the col names prove to be different across studies (especially filter pass)
        if column_map['filter_col'] in df.columns:
            df = df[df[column_map['filter_col']] == True]
        if column_map['alt_col'] in df.columns:
            df = df[df[column_map['alt_col']] != "<*>"]
        if "ALT_2" in df.columns and "ALT_3" in df.columns:
            df = df.drop(columns=["ALT_2", "ALT_3"])

        print(df.head())  # Debugging

        # Create unique variant identifier
        df["Variant"] = df["CHROM"].astype(str) + "." + df["POS"].astype(str) + "." + df["REF"] + "." + df["ALT"]

        # Add sample ID column
        df["Sample_ID"] = sample_id

        # Check if DataFrame is empty after filtering
        if df.empty:
            print(f"Warning: No variants left after filtering for sample {sample_id}")
            continue

        # Pivot to create binary presence-absence matrix
        df["Presence"] = 1  # Mark variants as present  
        variant_matrix = df.pivot(index="Sample_ID", columns="Variant", values="Presence").fillna(0)

        # Append to master dataframe
        master_df = pd.concat([master_df, variant_matrix], ignore_index=True)

    # Fill missing values with 0 and save to CSV
    master_df = master_df.fillna(0)
    master_df.to_csv(output_file, index=False)
    print(f"Successfully saved output to {output_file}")

# Parse command-line arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process VCF files into a binary presence-absence matrix.")
    parser.add_argument("--dir", required=True, help="Root directory containing SRX subdirectories")
    parser.add_argument("--output", required=True, help="Output CSV file path")
    
    # Arguments for column names
    parser.add_argument("--chrom", required=True, help="Column name for chromosome")
    parser.add_argument("--pos", required=True, help="Column name for position")
    parser.add_argument("--ref", required=True, help="Column name for reference allele")
    parser.add_argument("--alt", required=True, help="Column name for alternate allele")
    parser.add_argument("--filter", required=False, help="Column name for filter (default is 'FILTER_PASS')", default="FILTER_PASS")
    parser.add_argument("--samplesheet", required=True, help="Sample sheet containing the list of relevant SRX numbers")


    # Parse arguments
    args = parser.parse_args()

    # Create a column mapping dictionary
    column_mapping = {
        "chrom_col": args.chrom,
        "pos_col": args.pos,
        "ref_col": args.ref,
        "alt_col": args.alt,
        "filter_col": args.filter
    }
    
    # Run processing function
    process_vcf(args.dir, args.output, column_mapping)


'''
# Example Usage
python process_vcf.py \
  --dir "/ddn_exa/campbell/vpergola/Data/Tumour/wes/song2022/processed/variant_calling" \
  --output "variant_matrix.csv" \
  --chrom "CHROM" \
  --pos "POS" \
  --ref "REF" \
  --alt "ALT" \
  --filter "FILTER_PASS" \
  --samplesheet "" # TODO
'''

