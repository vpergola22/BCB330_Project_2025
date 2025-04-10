import pandas as pd
import argparse

def process_tsv(input_file, output_file, filter_file):
    # Read TSV file
    df = pd.read_csv(input_file, sep='\t')

    # Read filtered IDs (CSV file containing SRX_num)
    filtered_df = pd.read_csv(filter_file)

    # Extract the SRX_num from the filtered file
    filtered_srx = filtered_df['SRX_num'].tolist()

    # Filter the TSV based on SRX_num being in the filtered list
    df_filtered = df[df['experiment_accession'].isin(filtered_srx)]

    # Modify paths
    df_filtered['R1_Path'] = '/ddn_exa/campbell/vpergola/Data/Tumour/wes/song2022/fastq/' + df_filtered['id'].astype(str) + '_1.fastq.gz'
    df_filtered['R2_Path'] = '/ddn_exa/campbell/vpergola/Data/Tumour/wes/song2022/fastq/' + df_filtered['id'].astype(str) + '_2.fastq.gz'

    # Add a new column 'lane' with the same value and select relevant columns
    df_filtered['lane'] = 'lane_1'

    # Select relevant columns and save to CSV
    df_selected = df_filtered[['sample_accession', 'experiment_accession', 'lane', 'R1_Path', 'R2_Path']]
    
    # Rename the columns
    df_selected = df_selected.rename(columns={
        'sample_accession': 'patient',
        'experiment_accession': 'sample',
        'lane': 'lane',
        'R1_Path': 'fastq_1',
        'R2_Path': 'fastq_2'
    })

    # Save to CSV
    df_selected.to_csv(output_file, index=False)

    print(f"Sample sheet generated: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate sample sheet from TSV.")
    parser.add_argument("--input", required=True, help="Path to input TSV file")
    parser.add_argument("--output", required=True, help="Path to output CSV file")
    parser.add_argument("--filter", required=True, help="Path to filtered_id CSV file")

    args = parser.parse_args()

    process_tsv(args.input, args.output, args.filter)
