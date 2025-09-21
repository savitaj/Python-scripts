import scvi
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import pandas as pd

# Set up logging and figure parameters
sc.settings.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()

# --- 1. Data Loading and Preprocessing ---
print("1. Loading and preprocessing data...")

# Download and load the 10k PBMC dataset
# This is an excellent test dataset for single-cell analysis.
adata = scvi.data.pbmc_10k_protein_v3()

# Make a copy of the raw counts for later use
adata.layers["counts"] = adata.X.copy()

# Basic filtering and normalization
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata # Store the raw data for differential expression later

# --- 2. Simulating Batch Information (for demonstration) ---
# In a real-world scenario, you would have a 'batch' column
# in adata.obs (e.g., from different experiments or donors).
print("2. Simulating batch information for demonstration...")
adata.obs['batch'] = np.random.choice(['Batch_A', 'Batch_B', 'Batch_C'], size=adata.n_obs, p=[0.4, 0.4, 0.2])

# --- 3. Setup AnnData for scvi-tools (with Batch Correction) ---
print("3. Setting up AnnData object for scvi-tools with batch correction...")
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")

# --- 4. Model Training ---
print("4. Training the SCVI-based VAE model...")
# Initialize the SCVI model, a type of VAE tailored for single-cell data.
# n_latent: dimensionality of the latent space.
model = scvi.model.SCVI(adata, n_latent=20)

# Train the model. This is the core AI/ML step.
model.train(max_epochs=400)

# Save the trained model to a file
model.save("scvi_pbmc_model/", overwrite=True)

# --- 5. Latent Space Representation & Downstream Analysis ---
print("5. Performing downstream analysis: UMAP and Clustering...")

# Get the latent representation (the compressed, denoised data)
adata.obsm["X_scVI"] = model.get_latent_representation()

# Now, use this new latent representation for classic single-cell analysis tasks.
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata, min_dist=0.3)
sc.tl.leiden(adata, key_added="leiden_scvi")

# --- 6. Visualization ---
print("6. Generating visualizations...")

# Visualize the UMAP plot colored by Leiden clusters
fig, axs = plt.subplots(1, 2, figsize=(14, 6))
sc.pl.umap(adata, color="leiden_scvi", legend_loc="on data", title="UMAP with scVI Latent Space", ax=axs[0], show=False)
sc.pl.umap(adata, color="batch", title="UMAP Colored by Batch", ax=axs[1], show=False)
fig.suptitle('UMAP with SCVI Latent Space and Batch Correction', fontsize=16)
plt.tight_layout()
plt.show()

# Visualize the expression of key marker genes to annotate cell types.
marker_genes = ['MS4A1', 'CD3D', 'NKG7', 'PPBP', 'CD14']
sc.pl.umap(adata, color=marker_genes, cmap="viridis", title="Marker Gene Expression")
plt.show()

# --- 7. Differential Expression Analysis (with the SCVI model) ---
print("7. Performing differential expression analysis...")

# Find marker genes for the identified clusters.
# This method uses the scVI model's posterior distribution for a more robust test.
de_results = model.differential_expression(
    groupby="leiden_scvi",
    group1="0",  # Example: Find genes specific to cluster 0.
    mode="change"
)

# Display the top 10 differentially expressed genes for cluster 0.
print("Top 10 differentially expressed genes for cluster 0:")
print(de_results.sort_values("lfc_mean", ascending=False).head(10))

# --- 8. Gene Set Enrichment Analysis (GSEA) Input ---
# While not a full GSEA, this section prepares the data for an external tool.
print("\n8. Preparing data for Gene Set Enrichment Analysis...")

# Run Scanpy's differential expression to get ranked genes per cluster
sc.tl.rank_genes_groups(adata, 'leiden_scvi', method='wilcoxon')

# Get the top 100 genes for cluster 0
top_genes_cluster_0 = adata.uns['rank_genes_groups']['names']['0'][:100].tolist()

print(f"Top 100 genes for Cluster 0: {top_genes_cluster_0}")
print("\nThis gene list can be used with a GSEA tool (e.g., Enrichr) to find enriched pathways.")

print("\nAnalysis complete. The AnnData object `adata` contains all the results.")
print(adata)
