import pandas as pd
from io import StringIO


# .concat --> use to stack df into one row

def get_tcr_combs(df):
    # Filter to include only contig_1 and contig_2
    filtered_df = df[df['contig_id'].str.contains('contig_1|contig_2')]

    # Remove barcodes that appear more than twice
    barcode_counts = filtered_df['barcode'].value_counts()
    valid_barcodes = barcode_counts[barcode_counts <= 2].index
    filtered_df = filtered_df[filtered_df['barcode'].isin(valid_barcodes)]

    # Initialize a dictionary to store results and invalid results for each barcode
    result = {}
    invalid_barcodes = set()

    # Iterate through the remaining unique barcodes
    for barcode in filtered_df['barcode'].unique():
        subset = filtered_df[filtered_df['barcode'] == barcode]
        tcr_dict = {'TRA': {}, 'TRB': {}}
        
        tra_count = 0
        trb_count = 0
        is_valid = True  # Assume barcode is valid initially
        
        for _, row in subset.iterrows():
            if row['chain'] == 'TRA':
                if tra_count == 0:
                    tra_count += 1
                    tcr_dict['TRA']['V'] = row['v_gene']
                    tcr_dict['TRA']['J'] = row['j_gene']
                else:
                    # More than one TRA, mark as invalid
                    is_valid = False
                    break

            elif row['chain'] == 'TRB':
                if trb_count == 0:
                    trb_count += 1
                    tcr_dict['TRB']['V'] = row['v_gene']
                    tcr_dict['TRB']['J'] = row['j_gene']
                else:
                    # More than one TRB, mark as invalid
                    is_valid = False
                    break

        if is_valid:
            # Combine the results into the desired format
            result[barcode] = f"{tcr_dict['TRA'].get('V', 'NA')}|{tcr_dict['TRA'].get('J', 'NA')}|" \
                            f"{tcr_dict['TRB'].get('V', 'NA')}|{tcr_dict['TRB'].get('J', 'NA')}"
        else:
            # Add the barcode to the invalid list
            invalid_barcodes.add(barcode)

    # Filter out invalid barcodes
    filtered_result = {k: v for k, v in result.items() if k not in invalid_barcodes}

    # Convert the result to a DataFrame
    result_df = pd.DataFrame.from_dict(filtered_result, orient='index', columns=['TCR_combination'])
    result_df.reset_index(inplace=True)

    return result_df


def process_csv(input_file, GSM_num):
    # get data
    df = pd.read_csv(input_file)
    df.columns = df.columns.str.strip()

    # Get valid TCR combinations
    result_df = get_tcr_combs(df)

    # Remove rows with 'NA' or 'nan' values
    result_df = result_df[~result_df['TCR_combination'].str.contains('NA')]
    result_df = result_df[~result_df['TCR_combination'].str.contains('nan')]

    # Pivot to create binary presence-absence matrix
    binary_df = result_df.pivot_table(index='index', columns='TCR_combination', aggfunc='size', fill_value=0)

    # Reset index for saving
    binary_df.reset_index(inplace=True)

    # give all the barcodes the same value for GSM_num
    binary_df["GSM_num"] = GSM_num

    # get rid of barcodes
    binary_df = binary_df.drop('index', axis=1)

    # combine all the rows
    grouped = binary_df.groupby('GSM_num').max()
    grouped = grouped.reset_index()

    # Save the results
    return grouped
    # grouped.to_csv(output_file, index=False)



input_csv_list = [
    "PT11_Tumor_scVDJ/filtered_contig_annotations.csv",
    "PT56_Tumor_scVDJ/filtered_contig_annotations.csv",
    "PT55_Tumor_scVDJ/filtered_contig_annotations.csv",
    "PT53_Tumor_scVDJ/filtered_contig_annotations.csv",
    "PT52_Tumor_scVDJ/filtered_contig_annotations.csv",
    "PT50_Tumor_scVDJ/filtered_contig_annotations.csv",
    "PT47_Tumor_scVDJ/filtered_contig_annotations.csv",
    "PT35_Tumor_scVDJ/filtered_contig_annotations.csv"
]

# output_csv = "GSM8390840_binary_matrix.csv"

GSM_list = [
    "GSM8390840",
    "GSM8390826",
    "GSM8390828",
    "GSM8390830",
    "GSM8390832",
    "GSM8390834",
    "GSM8390836",
    "GSM8390838",
]


master_df = pd.DataFrame(columns=['GSM_num'])


# iterate through all the csv files from the study
for i in range(len(input_csv_list)):
    input_csv = "/ddn_exa/campbell/vpergola/Data/Tcell/sc/song2022/" + input_csv_list[i]
    GSM_num = GSM_list[i]
    
    curr = process_csv(input_csv, GSM_num)

    # add any new columns
    new_columns = {col: [0] * len(master_df) for col in curr.columns if col != 'GSM_num' and col not in master_df.columns}
    # Convert the dictionary of new columns to a DataFrame
    new_columns_df = pd.DataFrame(new_columns)
    # Concatenate the new columns with the existing master DataFrame
    master_df = pd.concat([master_df, new_columns_df], axis=1)
    # Add the new row to the master dataframe
    master_df = pd.concat([master_df, curr], ignore_index=True)

    # Update the TCR combination columns where the value is 1
    for col in curr.columns:
        if col != 'GSM_num':
            master_df.loc[master_df['GSM_num'] == curr['GSM_num'][0], col] = curr[col][0]


# Replace all blank values with 0
master_df = master_df.fillna(0)

# Convert the values to integers if required
# master_df = master_df.astype(int) TODO changing to ints caused an error

# Save the updated master dataframe to CSV
master_df.to_csv('master_data.csv', index=False)

