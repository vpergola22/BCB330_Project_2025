import allel as allel
import pandas as pd
import os
import argparse
import re

def process_vcf(directory, output_file, column_map, samplesheet):
    """
    Processes multiple VCF files in subdirectories and generates a binary presence-absence matrix.

    Parameters:
        directory (str): Directory where SRX subdirectories are located.
        output_file (str): Output CSV file path.
        column_map (dict): Dictionary containing column names for chrom, pos, ref, alt, and other necessary columns.
        samplesheet (str): Path to CSV file containing the list of SRX numbers to process.
    """
    # Load allowed SRX numbers from the samplesheet
    try:
        filtered_df = pd.read_csv(samplesheet)
        allowed_srx_numbers = set(filtered_df["SRX_num"].astype(str))  # Convert to a set for fast lookup
        print(f"Loaded {len(allowed_srx_numbers)} SRX numbers from {samplesheet}")
    except Exception as e:
        print(f"Error reading {samplesheet}: {e}")
        return

    # Initialize an empty list to store VCF file paths
    vcf_files = []

    # Recursively search for VCF files in SRX subdirectories
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith("deepvariant.vcf.gz"):  # Ensure correct VCF file format
                match = re.search(r"(SRX\d+)", root)  # Extract SRX number from folder path
                if match:
                    srx_num = match.group(1)
                    if srx_num in allowed_srx_numbers:  # Only add if in filtered list
                        vcf_files.append(os.path.join(root, file))

    print(f"Found {len(vcf_files)} VCF files matching the filtered SRX numbers.")

    master_df = pd.DataFrame()

    for vcf_file in vcf_files:
        # Extract SRX number again from filename (redundant but safe)
        match = re.search(r"(SRX\d+)", vcf_file)
        if match:
            srx_num = match.group(1)
        else:
            print(f"Warning: Could not extract SRX_num from {vcf_file}")
            continue  # Skip if no SRX match

        # Convert VCF to DataFrame
        try:
            df = allel.vcf_to_dataframe(vcf_file)
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

        # Apply filters (Modify if needed)
        if column_map['filter_col'] in df.columns:
            df = df[df[column_map['filter_col']] == True]
        if column_map['alt_col'] in df.columns:
            df = df[df[column_map['alt_col']] != "<*>"]
        if "ALT_2" in df.columns and "ALT_3" in df.columns:
            df = df.drop(columns=["ALT_2", "ALT_3"])

        # Create unique variant identifier
        df["Variant"] = df["CHROM"].astype(str) + "." + df["POS"].astype(str) + "." + df["REF"] + "." + df["ALT"]

        # Add SRX number as a column
        df["SRX_num"] = srx_num

        # Check if DataFrame is empty after filtering
        if df.empty:
            print(f"Warning: No variants left after filtering for sample {srx_num}")
            continue

        # Pivot to create binary presence-absence matrix
        df["Presence"] = 1  # Mark variants as present  
        variant_matrix = df.pivot(index="SRX_num", columns="Variant", values="Presence").fillna(0)

        # Append to master dataframe
        master_df = pd.concat([master_df, variant_matrix], ignore_index=False)

    # Fill missing values with 0 and save to CSV
    master_df = master_df.fillna(0)
    master_df.to_csv(output_file, index=True)
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
    process_vcf(args.dir, args.output, column_mapping, args.samplesheet)


"""
Example Usage:
python process_vcf.py \
  --dir "/ddn_exa/campbell/vpergola/Data/Tumour/wes/song2022/processed/variant_calling" \
  --output "variant_matrix.csv" \
  --chrom "CHROM" \
  --pos "POS" \
  --ref "REF" \
  --alt "ALT" \
  --filter "FILTER_PASS" \
  --samplesheet "filtered_ids.csv"
"""