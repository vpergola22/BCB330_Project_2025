import pandas as pd

# Load the MAF file into a DataFrame
maf_file_path = '/ddn_exa/campbell/vpergola/Data/Pipeline/meskit/final_mafs/SRX11758712.vep.maf'
maf_df = pd.read_csv(maf_file_path, sep='\t', comment='#')

# Extract the columns t_ref_count and t_alt_count
t_ref_alt_counts = maf_df[['VAF', 'Ref_allele_depth', 'Alt_allele_depth']]

# Display the first few rows to check the values
print(t_ref_alt_counts.head())
