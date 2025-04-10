import pandas as pd

def filter_tcr_data(input_file, output_file):
    """
    Filter TCR data from a CSV file based on specific criteria:
    1. Keep only contig_1 and contig_2
    2. Remove barcodes that appear more than twice
    3. Keep only barcodes with at most one TRA and one TRB
    
    Parameters:
    input_file (str): Path to the input CSV file
    output_file (str): Path to save the filtered CSV file
    """
    # Read the input CSV file
    df = pd.read_csv(input_file)
    df.columns = df.columns.str.strip()
    
    # Filter to include only contig_1 and contig_2
    filtered_df = df[df['contig_id'].str.contains('contig_1|contig_2')]
    
    # Remove barcodes that appear more than twice
    barcode_counts = filtered_df['barcode'].value_counts()
    valid_barcodes = barcode_counts[barcode_counts == 2].index
    filtered_df = filtered_df[filtered_df['barcode'].isin(valid_barcodes)]
    
    # Identify barcodes with invalid chain combinations
    invalid_barcodes = set()
    
    for barcode in filtered_df['barcode'].unique():
        subset = filtered_df[filtered_df['barcode'] == barcode]
        
        # Count TRA and TRB chains for this barcode
        tra_count = sum(subset['chain'] == 'TRA')
        trb_count = sum(subset['chain'] == 'TRB')
        
        # If more than one TRA or TRB, mark as invalid
        if tra_count > 1 or trb_count > 1:
            invalid_barcodes.add(barcode)
    
    # Remove invalid barcodes from the filtered dataframe
    final_df = filtered_df[~filtered_df['barcode'].isin(invalid_barcodes)]
    
    # Save the filtered data
    final_df.to_csv(output_file, index=False)
    
    return final_df

# Example usage
if __name__ == "__main__":
    # Process a single file
    # filter_tcr_data("input.csv", "filtered_output.csv")
    
    # Or process multiple files
    input_csv_list = [
        "PT11.csv",
        "PT56.csv",
        "PT55.csv",
        "PT53.csv",
        "PT52.csv",
        "PT50.csv",
        "PT47.csv",
        "PT35.csv"
    ]
    
    base_path = "/ddn_exa/campbell/vpergola/Data/Tcell/sc/song2022/immunarch/"
    
    for input_csv in input_csv_list:
        full_input_path = base_path + input_csv
        # Create output filename based on input
        sample_name = input_csv.split('/')[0]
        output_file = f"{sample_name}_filtered.csv"
        
        print(f"Processing {sample_name}...")
        filter_tcr_data(full_input_path, output_file)
        print(f"Saved filtered data to {output_file}")