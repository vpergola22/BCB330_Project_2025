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
mutation_data = pd.read_csv('/ddn_exa/campbell/vpergola/Data/Joined/song2022/song2022_maf.csv')
mutation_data.set_index('SRX_num', inplace=True)

# Load the ID mapping
id_mapping = pd.read_csv('/ddn_exa/campbell/vpergola/Data/Pipeline/ids/id_handling/filtered_ids.csv')

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

# Print bin sizes
for bin_name, bin_mutations in bins.items():
    print(f"Bin {bin_name}: {len(bin_mutations)} mutations")

# 3. Create results directories
os.makedirs('results', exist_ok=True)
os.makedirs('results/plots', exist_ok=True)

# Results will be stored here
all_results = []

# OPTIMIZATION: Fixed parameter combinations
# Removed None penalty to avoid warnings and separate l1/l2 penalties with different solvers
parameter_sets = [
    {'penalty': 'l1', 'C': 0.01, 'solver': 'liblinear'},
    {'penalty': 'l1', 'C': 0.1, 'solver': 'liblinear'},
    {'penalty': 'l1', 'C': 1, 'solver': 'liblinear'},
    {'penalty': 'l1', 'C': 10, 'solver': 'liblinear'},
    {'penalty': 'l2', 'C': 0.01, 'solver': 'liblinear'},
    {'penalty': 'l2', 'C': 0.1, 'solver': 'liblinear'},
    {'penalty': 'l2', 'C': 1, 'solver': 'liblinear'},
    {'penalty': 'l2', 'C': 10, 'solver': 'liblinear'},
    # No penalty version without C parameter to avoid warnings
    {'penalty': None, 'solver': 'newton-cg'}
]

# Function to calculate odds ratio and p-value from logistic regression model
def get_odds_ratio_and_pvalue(model, X, y):
    """Calculate odds ratios and p-values for logistic regression coefficients"""
    # Calculate odds ratios
    coef = model.coef_[0]
    odds_ratios = np.exp(coef)
    
    # Generate standard errors for coefficients using bootstrapping
    n_samples = X.shape[0]
    n_bootstraps = 1000
    bootstrapped_coefs = np.zeros((n_bootstraps, len(coef)))
    
    # Generate bootstrap samples and refit model
    for i in range(n_bootstraps):
        # Sample with replacement
        indices = np.random.choice(range(n_samples), size=n_samples, replace=True)
        X_boot, y_boot = X[indices], y[indices]
        
        try:
            # Refit model on bootstrap sample
            model_boot = LogisticRegression(
                penalty=model.penalty, 
                C=model.C if hasattr(model, 'C') else None,
                solver=model.solver,
                max_iter=1000,
                random_state=42
            )
            model_boot.fit(X_boot, y_boot)
            bootstrapped_coefs[i] = model_boot.coef_[0]
        except:
            # If fitting fails, reuse original coefficients
            bootstrapped_coefs[i] = coef
    
    # Calculate standard errors
    se = np.std(bootstrapped_coefs, axis=0)
    
    # Calculate z-scores and p-values
    z_scores = np.zeros_like(coef)
    p_values = np.ones_like(coef)
    
    # Only calculate where se > 0 to avoid division by zero
    valid_idx = se > 0
    z_scores[valid_idx] = coef[valid_idx] / se[valid_idx]
    p_values[valid_idx] = 2 * (1 - stats.norm.cdf(np.abs(z_scores[valid_idx])))
    
    return odds_ratios, p_values

# Function to process a single TCR
def process_tcr(tcr_idx, tcr_name, tcr_data, bins, mutation_data_pt, parameter_sets):
    results = []
    y = tcr_data[tcr_name].values
    
    # Skip TCRs with no variance
    if np.sum(y) == 0 or np.sum(y) == len(y):
        print(f"Skipping {tcr_name} - all samples have same value")
        return results
        
    print(f"Processing TCR {tcr_idx+1}: {tcr_name}")
    
    # Process each bin of mutations
    for bin_name, bin_mutations in bins.items():
        # Skip empty bins
        if not bin_mutations:
            continue
        
        # Get mutation data for this bin
        X = mutation_data_pt[bin_mutations].values
        
        # Skip if insufficient data or no variance
        if X.shape[0] < 2 or not np.any(np.std(X, axis=0) > 0):
            continue
        
        best_model = None
        best_score = -1
        best_params = {}
        best_metrics = {}
        
        # Try different model parameters
        for params in parameter_sets:
            try:
                # Set up model
                model = LogisticRegression(
                    penalty=params['penalty'],
                    C=params.get('C', 1.0),  # Default C=1.0 when not specified
                    solver=params['solver'],
                    max_iter=1000,
                    random_state=42
                )
                
                # Determine appropriate number of folds
                n_pos = np.sum(y)
                n_neg = len(y) - n_pos
                min_samples = min(n_pos, n_neg)
                
                # Check if we have enough samples for cross-validation
                if min_samples >= 2:
                    k = min(5, min_samples)
                    cv = StratifiedKFold(n_splits=k, shuffle=True, random_state=42)
                    
                    # Try cross-validation and fallback to train metrics if it fails
                    try:
                        cv_scores = cross_val_score(model, X, y, cv=cv, scoring='f1')
                        avg_score = np.mean(cv_scores)
                    except Exception as cv_err:
                        print(f"Cross-validation error, using training metrics: {str(cv_err)}")
                        model.fit(X, y)
                        y_pred = model.predict(X)
                        avg_score = f1_score(y, y_pred, zero_division=0)
                else:
                    # Not enough samples for cross-validation, use training metrics
                    model.fit(X, y)
                    y_pred = model.predict(X)
                    avg_score = f1_score(y, y_pred, zero_division=0)
                
                # Continue fitting and evaluation if not already done
                if not hasattr(model, 'coef_'):
                    model.fit(X, y)
                
                y_pred = model.predict(X)
                accuracy = accuracy_score(y, y_pred)
                f1 = f1_score(y, y_pred, zero_division=0)
                precision = precision_score(y, y_pred, zero_division=0)
                recall = recall_score(y, y_pred, zero_division=0)
                
                # Ensure we don't have perfect F1=1.0 for all cases (sign of overfitting)
                if f1 == 1.0:
                    # Check actual predictions vs ground truth
                    true_pos = np.sum(np.logical_and(y == 1, y_pred == 1))
                    if true_pos <= 1:  # Likely overfitting if only 1 or 0 true positives
                        continue
                
                if f1 > best_score:
                    # Calculate confusion matrix
                    cm = confusion_matrix(y, y_pred)
                    
                    # Get odds ratios and p-values
                    if np.any(model.coef_ != 0):
                        odds_ratios, p_values = get_odds_ratio_and_pvalue(model, X, y)
                        
                        best_model = model
                        best_score = f1
                        best_params = {
                            'penalty': params['penalty'],
                            'C': params.get('C', 'NA'),
                            'solver': params['solver']
                        }
                        best_metrics = {
                            'accuracy': accuracy,
                            'f1': f1,
                            'precision': precision,
                            'recall': recall,
                            'confusion_matrix': cm,
                            'odds_ratios': odds_ratios,
                            'p_values': p_values,
                            'cv_f1': avg_score
                        }
            
            except Exception as e:
                print(f"Error with parameters {params}: {str(e)}")
                continue
        
        # Record results if a good model was found
        if best_model is not None:
            # Collect results for significant features
            sig_features = []
            for i, (odds, p_val) in enumerate(zip(best_metrics['odds_ratios'], best_metrics['p_values'])):
                if p_val < 0.05:  # Significance threshold
                    sig_features.append({
                        'feature': bin_mutations[i],
                        'odds_ratio': odds,
                        'p_value': p_val,
                        'coefficient': best_model.coef_[0][i]
                    })
            
            result = {
                'tcr': tcr_name,
                'bin': bin_name,
                'params': best_params,
                'accuracy': best_metrics['accuracy'],
                'f1': best_metrics['f1'],
                'precision': best_metrics['precision'],
                'recall': best_metrics['recall'],
                'cv_f1': best_metrics['cv_f1'],
                'confusion_matrix': best_metrics['confusion_matrix'],
                'significant_features': sig_features,
                'num_mutations_in_bin': len(bin_mutations),
                'num_significant_features': len(sig_features)
            }
            
            results.append(result)
            
            # Plot confusion matrix
            plt.figure(figsize=(6, 5))
            sns.heatmap(best_metrics['confusion_matrix'], annot=True, fmt='d', cmap='Blues',
                       xticklabels=['Negative', 'Positive'], yticklabels=['Negative', 'Positive'])
            plt.title(f"Confusion Matrix - {tcr_name} - {bin_name}")
            plt.xlabel('Predicted')
            plt.ylabel('Actual')
            plt.tight_layout()
            plt.savefig(f"results/plots/cm_{tcr_name.replace('/', '_')}_{bin_name}.png")
            plt.close()
            
            # Save feature importance plot for significant features
            if len(sig_features) > 0:
                sig_features_sorted = sorted(sig_features, key=lambda x: abs(x['coefficient']), reverse=True)
                top_features = sig_features_sorted[:min(20, len(sig_features_sorted))]
                
                plt.figure(figsize=(10, 6))
                feature_names = [f['feature'] for f in top_features]
                odds_ratios = [f['odds_ratio'] for f in top_features]
                p_values = [f['p_value'] for f in top_features]
                
                # Create a colormap based on p-values
                colors = ['red' if p < 0.01 else 'orange' if p < 0.05 else 'gray' for p in p_values]
                
                plt.barh(range(len(feature_names)), odds_ratios, color=colors)
                plt.yticks(range(len(feature_names)), feature_names)
                plt.axvline(x=1, color='black', linestyle='--')
                plt.title(f"Significant Features - {tcr_name} - {bin_name}")
                plt.xlabel('Odds Ratio (log scale)')
                plt.xscale('log')
                plt.tight_layout()
                plt.savefig(f"results/plots/features_{tcr_name.replace('/', '_')}_{bin_name}.png")
                plt.close()
                
                # Save p-values table
                p_value_df = pd.DataFrame(sig_features_sorted)
                p_value_df.to_csv(f"results/plots/pvalues_{tcr_name.replace('/', '_')}_{bin_name}.csv", index=False)
    
    return results

# Main execution
# Process TCRs in parallel 
n_jobs = min(os.cpu_count(), 4)  # Use up to 4 cores
print(f"Processing {len(tcr_data.columns)} TCRs with {n_jobs} parallel jobs...")

# Process TCRs in batches to avoid memory issues
batch_size = 10
results_batches = []

for i in range(0, len(tcr_data.columns), batch_size):
    batch_end = min(i + batch_size, len(tcr_data.columns))
    tcr_batch = list(enumerate(tcr_data.columns[i:batch_end], i))
    
    print(f"Processing batch {i//batch_size + 1}/{(len(tcr_data.columns) + batch_size - 1)//batch_size}")
    
    # Process batch in parallel
    batch_results = Parallel(n_jobs=n_jobs)(
        delayed(process_tcr)(idx, name, tcr_data, bins, mutation_data_pt, parameter_sets) 
        for idx, name in tcr_batch
    )
    
    # Flatten results
    batch_results = [result for sublist in batch_results for result in sublist]
    results_batches.extend(batch_results)
    
    # Save intermediate results
    all_results.extend(batch_results)
    
    # Create intermediate dataframe and save
    if batch_results:
        results_df = pd.DataFrame([
            {
                'tcr': r['tcr'],
                'bin': r['bin'],
                'accuracy': r['accuracy'],
                'f1': r['f1'],
                'precision': r['precision'],
                'recall': r['recall'],
                'cv_f1': r['cv_f1'],
                'num_mutations_in_bin': r['num_mutations_in_bin'],
                'num_significant_features': r['num_significant_features'],
                'penalty': r['params']['penalty'],
                'C': r['params']['C'],
                'solver': r['params']['solver']
            }
            for r in all_results
        ])
        
        # Save intermediate results
        results_df.to_csv(f'results/intermediate_results_batch_{i//batch_size + 1}.csv', index=False)


# Add this code at the end of your script, replacing or enhancing your current visualization code

# Collate all results
if all_results:
    # Existing code for results_df creation...
    
    # Generate enhanced summary visualizations
    print("Generating enhanced summary visualizations...")
    
    # 1. Create a TCR summary dataframe with key metrics
    tcr_summary = []
    for tcr in results_df['tcr'].unique():
        tcr_data = results_df[results_df['tcr'] == tcr]
        
        summary = {
            'tcr_id': tcr,
            'mean_accuracy': tcr_data['accuracy'].mean(),
            'mean_f1': tcr_data['f1'].mean(),
            'mean_precision': tcr_data['precision'].mean(),
            'mean_recall': tcr_data['recall'].mean(),
            'mean_cv_f1': tcr_data['cv_f1'].mean(),
            'total_models': len(tcr_data),
            'significant_bins': len(tcr_data[tcr_data['num_significant_features'] > 0]),
            'best_bin': tcr_data.loc[tcr_data['f1'].idxmax(), 'bin'] if not tcr_data.empty else None,
            'best_f1': tcr_data['f1'].max() if not tcr_data.empty else 0
        }
        tcr_summary.append(summary)
    
    summary_df = pd.DataFrame(tcr_summary)
    summary_df.to_csv('results/tcr_summary.csv', index=False)
    
    # 2. Create comprehensive summary visualization
    plt.figure(figsize=(20, 16))
    
    # Plot 1: Compare TCRs by accuracy
    plt.subplot(3, 2, 1)
    if not summary_df.empty:
        summary_df_sorted = summary_df.sort_values("mean_accuracy", ascending=False)
        plt.bar(summary_df_sorted["tcr_id"].str.replace('/', '_'), summary_df_sorted["mean_accuracy"])
        plt.xticks(rotation=90)
        plt.xlabel("TCR ID")
        plt.ylabel("Mean Accuracy")
        plt.title("TCRs by Mean Accuracy")
        plt.grid(True, alpha=0.3)
    
    # Plot 2: Compare TCRs by significant bins
    plt.subplot(3, 2, 2)
    if not summary_df.empty:
        summary_df_sorted = summary_df.sort_values("significant_bins", ascending=False)
        plt.bar(summary_df_sorted["tcr_id"].str.replace('/', '_'), summary_df_sorted["significant_bins"])
        plt.xticks(rotation=90)
        plt.xlabel("TCR ID")
        plt.ylabel("Number of Significant Bins")
        plt.title("TCRs by Number of Significant Bins (with significant features)")
        plt.grid(True, alpha=0.3)
    
    # Plot 3: Heatmap of accuracy across all TCRs and bins
    plt.subplot(3, 2, 3)
    try:
        # Create the pivot table
        pivot_df = results_df.pivot(index="tcr", columns="bin", values="accuracy")
        
        # Create the heatmap
        sns.heatmap(pivot_df, cmap='viridis', annot=True, fmt=".2f", linewidths=.5)
        plt.xlabel("Mutation Bins")
        plt.ylabel("TCR ID")
        plt.title("Accuracy Heatmap: TCRs vs Mutation Bins")
    except Exception as e:
        print(f"Error creating accuracy heatmap: {str(e)}")
        plt.text(0.5, 0.5, "Insufficient data for accuracy heatmap", 
                 horizontalalignment='center', verticalalignment='center')
    
    # Plot 4: Heatmap of F1 scores across all TCRs and bins
    plt.subplot(3, 2, 4)
    try:
        # Create the pivot table for F1
        pivot_f1_df = results_df.pivot(index="tcr", columns="bin", values="f1")
        
        # Create the heatmap
        sns.heatmap(pivot_f1_df, cmap='plasma', annot=True, fmt=".2f", linewidths=.5)
        plt.xlabel("Mutation Bins")
        plt.ylabel("TCR ID")
        plt.title("F1 Score Heatmap: TCRs vs Mutation Bins")
    except Exception as e:
        print(f"Error creating F1 heatmap: {str(e)}")
        plt.text(0.5, 0.5, "Insufficient data for F1 heatmap", 
                 horizontalalignment='center', verticalalignment='center')
    
    # Plot 5: Heatmap of significant features across all TCRs and bins
    plt.subplot(3, 2, 5)
    try:
        # Create the pivot table for significant features
        pivot_sig_df = results_df.pivot(index="tcr", columns="bin", values="num_significant_features")
        
        # Create the heatmap
        sns.heatmap(pivot_sig_df, cmap='YlOrRd', annot=True, fmt="d", linewidths=.5)
        plt.xlabel("Mutation Bins")
        plt.ylabel("TCR ID")
        plt.title("Number of Significant Features: TCRs vs Mutation Bins")
    except Exception as e:
        print(f"Error creating significant features heatmap: {str(e)}")
        plt.text(0.5, 0.5, "Insufficient data for significant features heatmap", 
                 horizontalalignment='center', verticalalignment='center')
    
    # Plot 6: Box plot of F1 scores by bin
    plt.subplot(3, 2, 6)
    if not results_df.empty:
        sns.boxplot(x='bin', y='f1', data=results_df)
        plt.title('Distribution of F1 Scores by Mutation Bin')
        plt.ylabel('F1 Score')
        plt.xlabel('Mutation Bin')
        plt.xticks(rotation=45)
    
    plt.tight_layout()
    plt.savefig("results/tcr_summary_visualization.png", dpi=300)
    print("Enhanced summary visualization saved as tcr_summary_visualization.png")
    
    # Additional summary plots
    
    # 1. Bin comparison plot for each metric
    for metric in ['accuracy', 'f1', 'precision', 'recall', 'cv_f1']:
        plt.figure(figsize=(12, 8))
        sns.boxplot(x='bin', y=metric, data=results_df)
        plt.title(f'Distribution of {metric.upper()} by Mutation Bin')
        plt.ylabel(metric.upper())
        plt.xlabel('Mutation Bin')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(f"results/plots/bin_comparison_{metric}.png", dpi=300)
        plt.close()
    
    # 2. Top-performing TCRs
    top_tcrs = summary_df.nlargest(10, 'mean_f1')
    plt.figure(figsize=(12, 6))
    sns.barplot(x='tcr_id', y='mean_f1', data=top_tcrs)
    plt.title('Top 10 TCRs by Mean F1 Score')
    plt.ylabel('Mean F1 Score')
    plt.xlabel('TCR ID')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig("results/plots/top_tcrs.png", dpi=300)
    plt.close()
    
    # 3. Performance comparison across bins for specific TCRs
    top_5_tcrs = summary_df.nlargest(5, 'mean_f1')['tcr_id'].tolist()
    if top_5_tcrs:
        plt.figure(figsize=(15, 10))
        top_tcr_data = results_df[results_df['tcr'].isin(top_5_tcrs)]
        
        for i, metric in enumerate(['accuracy', 'f1', 'precision', 'recall']):
            plt.subplot(2, 2, i+1)
            sns.lineplot(x='bin', y=metric, hue='tcr', data=top_tcr_data, marker='o')
            plt.title(f'{metric.upper()} Across Bins for Top TCRs')
            plt.xticks(rotation=45)
            plt.legend(title='TCR ID')
        
        plt.tight_layout()
        plt.savefig("results/plots/top_tcrs_bin_comparison.png", dpi=300)
        plt.close()
    
    # 4. Correlation matrix of performance metrics
    plt.figure(figsize=(10, 8))
    correlation_metrics = ['accuracy', 'f1', 'precision', 'recall', 'cv_f1', 'num_significant_features']
    correlation_matrix = results_df[correlation_metrics].corr()
    sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', linewidths=.5)
    plt.title('Correlation Between Performance Metrics')
    plt.tight_layout()
    plt.savefig("results/plots/metric_correlations.png", dpi=300)
    plt.close()

# Show total runtime
end_time = time.time()
print(f"Analysis complete! Results saved to 'results' directory.")
print(f"Total runtime: {(end_time - start_time) / 60:.2f} minutes")