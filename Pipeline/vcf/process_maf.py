import pandas as pd
import os
import argparse
import re

def process_maf(directory, output_file, samplesheet):
    """
    Processes multiple MAF files and generates a binary presence-absence matrix based on the "dbSNP_RS" column.

    Parameters:
        directory (str): Directory where MAF files are located.
        output_file (str): Output CSV file path.
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

    # Initialize an empty dataframe
    master_df = pd.DataFrame()
    
    # Find and process MAF files
    maf_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".maf")]
    print(f"Found {len(maf_files)} MAF files.")
    
    for maf_file in maf_files:
        match = re.search(r"(SRX\d+)", maf_file)
        if match:
            srx_num = match.group(1)
        else:
            print(f"Warning: Could not extract SRX_num from {maf_file}")
            continue

        if srx_num not in allowed_srx_numbers:
            continue

        try:
            df = pd.read_csv(maf_file, sep='\t', header=1, low_memory=False)
        except Exception as e:
            print(f"Error processing {maf_file}: {e}")
            continue
        
        print(df.columns)

        if "dbSNP_RS" not in df.columns:
            print(f"Skipping {maf_file}: 'dbSNP_RS' column not found.")
            continue
        
        df = df.dropna(subset=["dbSNP_RS"])  # Remove rows with missing RSIDs
        df["dbSNP_RS"] = df["dbSNP_RS"].astype(str)  # Ensure RSIDs are strings
        
        df["Presence"] = 1
        variant_matrix = df.pivot_table(index="dbSNP_RS", columns=[], values="Presence", aggfunc='max').T
        variant_matrix["SRX_num"] = srx_num
        
        master_df = pd.concat([master_df, variant_matrix], ignore_index=False)
    
    master_df = master_df.fillna(0)
    master_df.set_index("SRX_num", inplace=True)
    master_df.to_csv(output_file)
    print(f"Successfully saved output to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process MAF files into a binary presence-absence matrix using dbSNP_RS.")
    parser.add_argument("--dir", required=True, help="Directory containing MAF files")
    parser.add_argument("--output", required=True, help="Output CSV file path")
    parser.add_argument("--samplesheet", required=True, help="Sample sheet containing relevant SRX numbers")
    
    args = parser.parse_args()
    process_maf(args.dir, args.output, args.samplesheet)
