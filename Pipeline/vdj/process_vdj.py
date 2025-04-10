import pandas as pd
import os
import argparse
import re


def find_vdj_files(base_dir):
    """
    Finds GSM numbers from ZIP filenames and corresponding CSV paths.

    Parameters:
        base_dir (str): Directory containing study data.

    Returns:
        list: List of paths to relevant CSV files.
        list: List of corresponding GSM numbers.
    """
    gsm_mapping = {}

    # Extract GSM numbers from ZIP filenames
    for filename in os.listdir(base_dir):
        match = re.match(r"(GSM\d+)_PT\d+_Tumor_scVDJ\.zip", filename)
        if match:
            gsm_number = match.group(1)
            patient_id = filename.split("_")[1]  # Extract PTXX
            gsm_mapping[patient_id] = gsm_number

    # Find corresponding CSV files in subdirectories
    input_files = []
    gsm_list = []

    for subdir in os.listdir(base_dir):
        subdir_path = os.path.join(base_dir, subdir)

        if os.path.isdir(subdir_path) and "Tumor_scVDJ" in subdir:
            csv_file = os.path.join(subdir_path, "filtered_contig_annotations.csv")

            if os.path.exists(csv_file):
                patient_id = subdir.split("_")[0]  # Extract PTXX
                if patient_id in gsm_mapping:
                    input_files.append(csv_file)
                    gsm_list.append(gsm_mapping[patient_id])

    return input_files, gsm_list


def get_tcr_combs(df, chain_col, v_col, j_col, barcode_col, contig_col):
    """
    Extracts unique TCR combinations for each barcode while filtering out invalid ones.
    """
    filtered_df = df[df[contig_col].str.contains('contig_1|contig_2', na=False)]
    barcode_counts = filtered_df[barcode_col].value_counts()
    valid_barcodes = barcode_counts[barcode_counts <= 2].index
    filtered_df = filtered_df[filtered_df[barcode_col].isin(valid_barcodes)]

    result = {}
    invalid_barcodes = set()

    for barcode in filtered_df[barcode_col].unique():
        subset = filtered_df[filtered_df[barcode_col] == barcode]
        tcr_dict = {'TRA': {}, 'TRB': {}}
        
        tra_count = 0
        trb_count = 0
        is_valid = True  

        for _, row in subset.iterrows():
            chain = row[chain_col]
            v_gene = row[v_col]
            j_gene = row[j_col]

            if chain == 'TRA':
                if tra_count == 0:
                    tcr_dict['TRA']['V'] = v_gene
                    tcr_dict['TRA']['J'] = j_gene
                    tra_count += 1
                else:
                    is_valid = False
                    break

            elif chain == 'TRB':
                if trb_count == 0:
                    tcr_dict['TRB']['V'] = v_gene
                    tcr_dict['TRB']['J'] = j_gene
                    trb_count += 1
                else:
                    is_valid = False
                    break

        if is_valid:
            result[barcode] = f"{tcr_dict['TRA'].get('V', 'NA')}.{tcr_dict['TRA'].get('J', 'NA')}." \
                              f"{tcr_dict['TRB'].get('V', 'NA')}.{tcr_dict['TRB'].get('J', 'NA')}"
        else:
            invalid_barcodes.add(barcode)

    filtered_result = {k: v for k, v in result.items() if k not in invalid_barcodes}
    result_df = pd.DataFrame.from_dict(filtered_result, orient='index', columns=['TCR_combination'])
    result_df.reset_index(inplace=True)
    result_df.rename(columns={'index': barcode_col}, inplace=True)

    return result_df


def process_csv(input_file, gsm_num, chain_col, v_col, j_col, barcode_col, contig_col):
    """
    Processes a CSV file and generates a binary presence-absence matrix.
    """
    df = pd.read_csv(input_file)
    df.columns = df.columns.str.strip()

    result_df = get_tcr_combs(df, chain_col, v_col, j_col, barcode_col, contig_col)
    result_df = result_df[~result_df['TCR_combination'].str.contains('NA.nan', na=False)]

    binary_df = result_df.pivot_table(index=barcode_col, columns='TCR_combination', aggfunc='size', fill_value=0)
    binary_df.reset_index(inplace=True)

    binary_df["GSM_num"] = gsm_num
    binary_df.drop(barcode_col, axis=1, inplace=True)

    grouped = binary_df.groupby('GSM_num').max().reset_index()

    return grouped


def process_dataset(base_dir, output_csv, chain_col, v_col, j_col, barcode_col, contig_col):
    """
    Automatically processes the dataset using detected files and GSM numbers.
    """
    input_files, gsm_list = find_vdj_files(base_dir)

    if not input_files:
        print("No valid files found.")
        return

    master_df = pd.DataFrame(columns=['GSM_num'])

    for i, input_csv in enumerate(input_files):
        gsm_num = gsm_list[i]
        print(f"Processing {input_csv} for {gsm_num}...")

        curr = process_csv(input_csv, gsm_num, chain_col, v_col, j_col, barcode_col, contig_col)

        new_columns = {col: [0] * len(master_df) for col in curr.columns if col != 'GSM_num' and col not in master_df.columns}
        master_df = pd.concat([master_df, pd.DataFrame(new_columns)], axis=1)
        master_df = pd.concat([master_df, curr], ignore_index=True)

        for col in curr.columns:
            if col != 'GSM_num':
                master_df.loc[master_df['GSM_num'] == gsm_num, col] = curr[col][0]

    master_df.fillna(0, inplace=True)
    master_df.to_csv(output_csv, index=False)
    print(f"Saved processed data to {output_csv}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process VDJ sequencing data.")
    parser.add_argument("base_dir", type=str, help="Path to study directory containing VDJ data.")
    parser.add_argument("output_csv", type=str, help="Output file path for the processed data.")
    parser.add_argument("--chain_col", type=str, default="chain", help="Column name for TCR chain.")
    parser.add_argument("--v_col", type=str, default="v_gene", help="Column name for V gene.")
    parser.add_argument("--j_col", type=str, default="j_gene", help="Column name for J gene.")
    parser.add_argument("--barcode_col", type=str, default="barcode", help="Column name for cell barcode.")
    parser.add_argument("--contig_col", type=str, default="contig_id", help="Column name for contig ID.")

    args = parser.parse_args()

    process_dataset(
        args.base_dir, args.output_csv, args.chain_col, args.v_col, args.j_col, args.barcode_col, args.contig_col
    )



'''
# Example Usage
python process_vdj.py /ddn_exa/campbell/vpergola/Data/Tcell/sc/song2022 /path/to/output_folder/processed_data.csv \
  --chain_col "chain" \
  --v_col "v_gene" \
  --j_col "j_gene" \
  --barcode_col "barcode" \
  --contig_col "contig_id"
'''

