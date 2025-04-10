import pandas as pd

def transform_tcr_data(input_csv, output_csv):
    # Load the CSV file
    df = pd.read_csv(input_csv)
    
    # Melt the DataFrame to long format
    df_melted = df.melt(id_vars=["V.name", "J.name"], var_name="Sample", value_name="Count")
    
    # Drop "Samples" row since it is not an actual sample
    df_melted = df_melted[df_melted["Sample"] != "Samples"]
    
    # Create new column names by combining V.name and J.name
    df_melted["VJ_combined"] = df_melted["V.name"] + "_" + df_melted["J.name"]
    
    # Pivot table to get one-hot encoding (presence = 1, absence = 0)
    df_pivot = df_melted.pivot_table(index="Sample", columns="VJ_combined", values="Count", aggfunc=lambda x: 1 if any(x.notna()) else 0)
    
    # Fill NA values with 0
    df_pivot = df_pivot.fillna(0).astype(int)
    
    # Reset index to have Sample as a column
    df_pivot.reset_index(inplace=True)
    
    # Save the transformed DataFrame to a new CSV file
    df_pivot.to_csv(output_csv, index=False)
    print(f"Transformed data saved to {output_csv}")

# Example usage
transform_tcr_data("public_tcrs_vj.csv", "public_tcrs_vj_new.csv")
transform_tcr_data("public_tcrs_vjaa.csv", "public_tcrs_vjaa_new.csv")
