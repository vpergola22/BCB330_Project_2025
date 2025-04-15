import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix, precision_score, recall_score
from sklearn.model_selection import cross_val_score, StratifiedKFold
from scipy import stats
import os
import time
from joblib import Parallel, delayed

# Start timing
start_time = time.time()

# 1. Data Loading and Preprocessing
print("Loading data...")

# Load the TCR data
tcr_data = pd.read_csv('public_tcrs_vj.csv')
tcr_data.set_index('Sample', inplace=True)

# Load the mutation data
mutation_data = pd.read_csv('/ddn_exa/campbell/vpergola/BCB330_Project_2025/Data/Joined/song2022/song2022_maf.csv')
mutation_data.set_index('SRX_num', inplace=True)

# Load the ID mapping
id_mapping = pd.read_csv('/ddn_exa/campbell/vpergola/BCB330_Project_2025/Pipeline/ids/id_handling/filtered_ids.csv')

# Create mappings more efficiently
srx_to_pt = dict(zip(id_mapping['SRX_num'], id_mapping['PT_num']))
pt_to_srx = {v: k for k, v in srx_to_pt.items()}  # Reverse mapping

# Filter mutation data to only include samples that are in the mapping
valid_srx = list(set(mutation_data.index) & set(srx_to_pt.keys()))
mutation_data = mutation_data.loc[valid_srx]

# Create mutation data with PT sample IDs
mutation_data_dict = {}
for srx in valid_srx:
    pt = srx_to_pt[srx]
    mutation_data_dict[pt] = mutation_data.loc[srx].to_dict()

# Convert dictionary to DataFrame at once
mutation_data_pt = pd.DataFrame.from_dict(mutation_data_dict, orient='index')

# Filter tcr_data to include only samples in the mapping
tcr_data = tcr_data[tcr_data.index.isin(pt_to_srx)]

# Get common samples and align dataframes
common_samples = sorted(list(set(tcr_data.index) & set(mutation_data_pt.index)))
tcr_data = tcr_data.loc[common_samples]
mutation_data_pt = mutation_data_pt.loc[common_samples]

print(f"Common samples: {len(common_samples)}")
print(f"TCR data shape: {tcr_data.shape}")
print(f"Mutation data shape: {mutation_data_pt.shape}")

# 2. Create mutation bins based on frequency
print("Creating mutation bins...")

# Count occurrences of each mutation across samples
mutation_counts = mutation_data_pt.sum()
print(f"Total mutations: {len(mutation_counts)}")

# Create bins
bin_exactly_1 = mutation_counts[mutation_counts == 1].index.tolist()
bin_exactly_8 = mutation_counts[mutation_counts == 8].index.tolist()

bins = {
    'exactly_1': bin_exactly_1,
    'exactly_8': bin_exactly_8
}

# Add cumulative bins from 2 to 7
for i in range(2, 8):
    bins[f'1to{i}'] = mutation_counts[(mutation_counts > 1) & (mutation_counts <= i)].index.tolist()

# Prepare data for plotting
bin_names = []
bin_sizes = []

for bin_name, bin_mutations in bins.items():
    bin_names.append(bin_name)
    bin_sizes.append(len(bin_mutations))

# Plotting
plt.figure(figsize=(10, 6))
sns.barplot(x=bin_names, y=bin_sizes)
plt.xlabel("Mutation Bin")
plt.ylabel("Number of Mutations")
plt.title("Number of Mutations per Bin")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("mutation_bins_barplot.png")  # Save the plot if desired
plt.show()
