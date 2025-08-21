import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# =====================
# Synthetic data maker
# =====================

"""
    This function SIMULATES gene expression data for a given species.
    
     IMPORTANT: This is NOT real biological data.
    It is a simulation that mimics what the actual tool would do with real input.
    
    - In Project 1, our tool cleaned real gene expression data and created 
      a co-expression matrix across species.
    - In Project 2, we extended the tool into a simulator that shows how 
      conditions of a gene could be visualized (heatmaps, patterns, etc.).
    """
def generate_synthetic_data(species, output_dir="data", num_genes=1000):
    os.makedirs(output_dir, exist_ok=True)

    filepath = os.path.join(output_dir, f"{species}.csv")
    if os.path.exists(filepath):
        print(f"[INFO] Found existing {filepath}, skipping synthetic generation.")
        return filepath

    np.random.seed(42)  # reproducibility
    genes = [f"Gene_{i}" for i in range(num_genes)]
    log2fc = np.random.normal(0, 2, num_genes)  # log2 fold changes
    pvalues = np.random.uniform(0.0001, 1, num_genes)  # p-values
    pvalues[np.random.choice(num_genes, size=50, replace=False)] = np.random.uniform(0.0001, 0.05, 50)  

    df = pd.DataFrame({
        "Gene": genes,
        "log2FoldChange": log2fc,
        "pvalue": pvalues
    })

    df.to_csv(filepath, index=False)
    print(f"[INFO] Synthetic data saved to {filepath}")
    return filepath


# =====================
# Volcano plot function
# =====================
def volcano_plot(species, data_file):
    df = pd.read_csv(data_file)

    df['-log10pvalue'] = -np.log10(df['pvalue'])

    plt.figure(figsize=(8,6))
    plt.scatter(df['log2FoldChange'], df['-log10pvalue'],
                c="grey", alpha=0.6, edgecolor='none')

    sig = (df['pvalue'] < 0.05) & (abs(df['log2FoldChange']) > 1)
    plt.scatter(df.loc[sig, 'log2FoldChange'], df.loc[sig, '-log10pvalue'],
                c="red", alpha=0.8, edgecolor='none')

    plt.axhline(-np.log10(0.05), linestyle='--', color='blue')
    plt.axvline(1, linestyle='--', color='green')
    plt.axvline(-1, linestyle='--', color='green')

    plt.title(f"Volcano Plot - {species}")
    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10(p-value)")
    plt.savefig(f"{species}_volcano.png", dpi=300)
    plt.close()
    print(f"[INFO] Volcano plot saved as {species}_volcano.png")


# =====================
# Main driver
# =====================
if __name__ == "__main__":
    species_list = input("Enter species names (comma-separated): ").split(",")

    for species in [s.strip() for s in species_list]:
        csv_path = generate_synthetic_data(species)
        volcano_plot(species, csv_path)
