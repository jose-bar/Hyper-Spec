import pandas as pd
import numpy as np
from scipy.stats import entropy
import sys
import os

# Check that both cluster size and runtime arguments were provided
if len(sys.argv) < 3:
    print("Usage: python results.py <cluster_size> <runtime>")
    sys.exit(1)

# Set cluster size and runtime from command line arguments
cluster_size = float(sys.argv[1])
runtime = float(sys.argv[2])

# Load the data
df = pd.read_parquet('output.csv.parquet')

# Define Metrics

def clustered_spectra_ratio(df):
    """Ratio of clustered spectra to total spectra."""
    clustered_count = df['cluster'].nunique()  # Number of unique clusters
    total_spectra = len(df)
    return clustered_count / total_spectra

def incorrect_clustering_ratio(df):
    """Ratio of incorrectly clustered spectra based on peptide labels."""
    incorrect_count = 0
    total_clustered_identified = 0

    for cluster_id, group in df.groupby('cluster'):
        # Get the most frequent peptide label in this cluster
        if group['identifier'].nunique() > 1:
            most_frequent_label = group['identifier'].mode()[0]
            # Count spectra with incorrect labels in this cluster
            incorrect_count += len(group[group['identifier'] != most_frequent_label])

        # Count total identified spectra in the cluster
        total_clustered_identified += len(group)

    return incorrect_count / total_clustered_identified if total_clustered_identified > 0 else 0

def completeness(df):
    """Compute the completeness based on entropy of cluster distribution for peptide assignments."""
    peptide_labels = df['identifier'].unique()

    # Entropy calculation for completeness
    H_C_given_P = 0  # Conditional entropy of cluster distribution given peptide assignments
    H_P = 0          # Entropy of peptide assignments

    for peptide_label in peptide_labels:
        peptide_group = df[df['identifier'] == peptide_label]
        cluster_counts = peptide_group['cluster'].value_counts(normalize=True)
        
        # Compute H_C_given_P for this peptide label
        H_C_given_P += entropy(cluster_counts) * (len(peptide_group) / len(df))
        
        # Contribution to H_P
        H_P += - (len(peptide_group) / len(df)) * np.log(len(peptide_group) / len(df))

    completeness_score = 1 - (H_C_given_P / H_P) if H_P > 0 else 1
    return completeness_score

# Calculate Metrics
metrics = {
    "Cluster Size": cluster_size,
    "Clustered Spectra Ratio": clustered_spectra_ratio(df),
    "Incorrect Clustering Ratio": incorrect_clustering_ratio(df),
    "Completeness": completeness(df),
    "Runtime": runtime
}

# File path to store the results
output_file = 'clustering_metrics.csv'

# Check if file exists to determine if we need to write headers
file_exists = os.path.isfile(output_file)

# Append the results to the CSV file
with open(output_file, 'a') as f:
    if not file_exists:
        # Write the header if file does not exist
        f.write("Cluster Size,Clustered Spectra Ratio,Incorrect Clustering Ratio,Completeness,Runtime\n")
    
    # Append the new results
    f.write(f"{metrics['Cluster Size']},{metrics['Clustered Spectra Ratio']:.4f},{metrics['Incorrect Clustering Ratio']:.4f},{metrics['Completeness']:.4f},{metrics['Runtime']:.4f}\n")

# Output Results
for metric, value in metrics.items():
    print(f"{metric}: {value:.4f}" if isinstance(value, float) else f"{metric}: {value}")

