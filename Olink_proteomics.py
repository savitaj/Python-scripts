
## Libraries to be Installed
## pip install pandas numpy scipy statsmodels seaborn matplotlib

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats import multitest
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Set up figure parameters for clean plots
sns.set_theme(style="whitegrid")
plt.rcParams['figure.dpi'] = 150
plt.rcParams['figure.figsize'] = (10, 8)

# --- 1. Data Loading and Preparation (with a Mock Dataset) ---
print("1. Loading and preparing data...")

# IMPORTANT: In a real analysis, replace this section with your data loading logic.
# Your Olink NPX data file (e.g., a CSV or Excel file) will contain NPX values
# and metadata columns (e.g., SampleID, PatientID, Group, etc.).

# Mock data generation for demonstration (replace with your own .csv or .xlsx Olink NPX file)
def generate_mock_olink_data():
    np.random.seed(42)
    
    # Create mock metadata
    groups = ['Control'] * 20 + ['Treated'] * 20
    sample_ids = [f'Sample_{i+1}' for i in range(40)]
    
    # Create mock NPX data for 100 proteins
    protein_ids = [f'Protein_{i+1}' for i in range(100)]
    
    # Generate random NPX values (log2 scale)
    npx_data = pd.DataFrame(np.random.normal(loc=np.tile([5, 5], 20), scale=0.8, size=(40, 100)),
                            index=sample_ids, columns=protein_ids)
    
    # Introduce a mock "treatment effect" for a few proteins
    npx_data['Protein_10'] = np.random.normal(loc=10, scale=0.8, size=40)
    npx_data['Protein_25'] = np.random.normal(loc=3, scale=0.8, size=40)
    npx_data['Protein_50'] = np.random.normal(loc=12, scale=0.8, size=40)
    
    # Introduce a mock "batch effect" to demonstrate correction
    batch_a = npx_data.iloc[:20, :] + 1 # Add 1 to NPX for batch A
    batch_b = npx_data.iloc[20:, :]
    npx_data = pd.concat([batch_a, batch_b])

    # Combine into a long format for easier analysis
    npx_long = npx_data.melt(ignore_index=False, var_name='Protein', value_name='NPX')
    npx_long.index.name = 'SampleID'
    npx_long = npx_long.reset_index()

    # Add metadata to the long-format data
    npx_long['Group'] = npx_long['SampleID'].apply(lambda x: 'Control' if int(x.split('_')[1]) <= 20 else 'Treated')
    npx_long['Batch'] = npx_long['SampleID'].apply(lambda x: 'Batch_A' if int(x.split('_')[1]) <= 20 else 'Batch_B')

    return npx_long

# Load the mock data
df_npx = generate_mock_olink_data()
print("Data loaded. Head of the dataframe:")
print(df_npx.head())

# --- 2. Quality Control and Visualization ---
print("\n2. Performing Quality Control and Visualization...")

# Plotting NPX distribution across all samples
plt.figure(figsize=(12, 6))
sns.boxplot(x='SampleID', y='NPX', data=df_npx)
plt.xticks(rotation=90)
plt.title('NPX Distribution Across Samples')
plt.show()

# Plotting NPX distribution across groups
plt.figure(figsize=(8, 6))
sns.boxplot(x='Group', y='NPX', data=df_npx)
plt.title('NPX Distribution Across Groups')
plt.show()

# Plotting NPX distribution across batches
plt.figure(figsize=(8, 6))
sns.boxplot(x='Batch', y='NPX', data=df_npx)
plt.title('NPX Distribution Across Batches')
plt.show()

# --- 3. Differential Expression Analysis (T-Test) ---
print("\n3. Performing Differential Expression Analysis...")

# Create a wide format dataframe for statistical testing
df_npx_wide = df_npx.pivot_table(index='SampleID', columns='Protein', values='NPX')
df_npx_wide['Group'] = df_npx['Group'].unique()
df_npx_wide['Batch'] = df_npx['Batch'].unique()
df_npx_wide = df_npx_wide.set_index([df_npx_wide.index, 'Group', 'Batch'])

# Split data into groups for comparison
group1_data = df_npx_wide.loc[df_npx_wide.index.get_level_values('Group') == 'Control']
group2_data = df_npx_wide.loc[df_npx_wide.index.get_level_values('Group') == 'Treated']

# Prepare a DataFrame to store results
results = pd.DataFrame(columns=['Protein', 'p_value', 'log2FC'])

# Perform Welch's t-test for each protein
for protein in df_npx_wide.columns:
    t_stat, p_val = stats.ttest_ind(
        group2_data[protein],
        group1_data[protein],
        equal_var=False  # Welch's t-test, more robust for unequal variances
    )
    
    # Calculate log2 fold change
    log2FC = np.mean(group2_data[protein]) - np.mean(group1_data[protein])
    
    # Append to results
    results.loc[len(results)] = [protein, p_val, log2FC]

# Correct for multiple testing using Benjamini-Hochberg (FDR)
results['FDR_p_value'] = multitest.multipletests(results['p_value'], method='fdr_bh')[1]

# Sort results by adjusted p-value
results = results.sort_values(by='FDR_p_value')
print("Differential expression analysis complete. Top 10 results:")
print(results.head(10))

# --- 4. Volcano Plot Visualization ---
print("\n4. Generating Volcano Plot...")

# Define significance threshold
alpha = 0.05
results['significant'] = results['FDR_p_value'] < alpha
results['neg_log10_p_value'] = -np.log10(results['FDR_p_value'])

plt.figure(figsize=(10, 10))
sns.scatterplot(
    x='log2FC',
    y='neg_log10_p_value',
    data=results,
    hue='significant',
    palette=['gray', 'red'],
    alpha=0.7
)

# Add labels for top proteins
for i, row in results.head(5).iterrows():
    plt.text(
        row['log2FC'] + 0.05,
        row['neg_log10_p_value'],
        row['Protein'],
        fontsize=9,
        ha='left',
        va='bottom'
    )

plt.axhline(y=-np.log10(alpha), color='r', linestyle='--', linewidth=1, label=f'FDR < {alpha}')
plt.title('Volcano Plot of Differential Protein Expression')
plt.xlabel('Log2 Fold Change (Treated vs. Control)')
plt.ylabel('-log10(FDR Adjusted p-value)')
plt.legend(title='Significant')
plt.show()

# --- 5. Heatmap Visualization of Top Proteins ---
print("\n5. Generating Heatmap...")

# Get the top N significant proteins
top_proteins = results[results['significant']].head(10)['Protein'].tolist()

if top_proteins:
    # Get NPX data for top proteins
    heatmap_data = df_npx.pivot_table(index='SampleID', columns='Protein', values='NPX')
    heatmap_data = heatmap_data[top_proteins]
    
    # Create annotations for the heatmap
    row_colors = df_npx.drop_duplicates(subset=['SampleID']).set_index('SampleID')[['Group', 'Batch']]
    group_color_map = {'Control': 'blue', 'Treated': 'orange'}
    batch_color_map = {'Batch_A': 'green', 'Batch_B': 'purple'}
    
    row_colors['Group'] = row_colors['Group'].map(group_color_map)
    row_colors['Batch'] = row_colors['Batch'].map(batch_color_map)

    # Plot the heatmap
    sns.clustermap(
        heatmap_data,
        row_colors=row_colors,
        cmap='viridis',
        figsize=(12, 12),
        dendrogram_ratio=(0.2, 0.2),
        cbar_pos=(0.02, 0.8, 0.05, 0.18)
    )
    plt.suptitle("Heatmap of Top Differentially Expressed Proteins", y=1.02)
    plt.show()

else:
    print("No significant proteins found to plot a heatmap.")

print("\nPipeline complete. The analysis results are stored in the 'results' DataFrame.")
