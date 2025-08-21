
import os
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

def clean_abundance_file(filepath, output_dir):
    """
    Load abundance.tsv, keep rows where est_counts > 0, save cleaned CSV.
    Returns cleaned dataframe and sample name.
    """
    try:
        df = pd.read_csv(filepath, sep='\t')
        df = df[df['est_counts'] > 0]
        sample = os.path.splitext(os.path.basename(filepath))[0]
        outpath = os.path.join(output_dir, f"{sample}_cleaned.csv")
        df.to_csv(outpath, index=False)
        return df, sample
    except Exception as e:
        print(f"[WARN] Skipping {filepath}: {e}")
        return None, None

def merge_samples(cleaned_files):
    """
    Merge cleaned abundance files into a single matrix (genes x samples).
    """
    dfs = []
    sample_names = []
    for f in cleaned_files:
        try:
            df = pd.read_csv(f)
            sample = os.path.basename(f).replace('_cleaned.csv', '')
            # Rename sample to rose or rice if possible
            if 'rose' in sample.lower():
                sample = 'abundance_rose'
            elif 'rice' in sample.lower():
                sample = 'abundance_rice'
            df = df[['target_id', 'est_counts']].copy()
            df = df.rename(columns={'est_counts': sample})
            df = df.set_index('target_id')
            dfs.append(df)
            sample_names.append(sample)
        except Exception as e:
            print(f"[WARN] Could not merge {f}: {e}")
    if not dfs:
        raise ValueError("No cleaned files to merge.")
    merged = pd.concat(dfs, axis=1).fillna(0)
    merged.to_csv('merged_matrix.csv')
    return merged

def visualize_heatmap(df, output_dir, top_n=50):
    """
    Plot heatmap of top variable genes (by variance).
    """
    var_genes = df.var(axis=1).sort_values(ascending=False).head(top_n).index
    data = df.loc[var_genes]
    # Create a more detailed heatmap
    plt.figure(figsize=(14, 10))
    ax = sns.heatmap(data, cmap='viridis', annot=True, fmt='.1f', linewidths=0.5, linecolor='gray', cbar_kws={'label': 'Abundance'})
    ax.set_xlabel('Sample')
    ax.set_ylabel('Gene')
    ax.set_title(f'Heatmap of Top {top_n} Variable Genes (Rose vs Rice)')
    plt.xticks(rotation=45)
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'detailed_heatmap_rose_vs_rice.png'))
    plt.close()
    print(f"Saved detailed heatmap of top {top_n} variable genes as detailed_heatmap_rose_vs_rice.png.")

# Visualize PCA of samples  
def visualize_pca(df, output_dir):

    # PCA plot of samples based on gene expression.

    pca = PCA(n_components=2)
    X = df.T.values
    pcs = pca.fit_transform(X)
    plt.figure(figsize=(8,6))
    plt.scatter(pcs[:,0], pcs[:,1], c='blue')
    for i, sample in enumerate(df.columns):
        plt.text(pcs[i,0], pcs[i,1], sample)
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('PCA of Samples')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'pca_plot.png'))
    plt.close()
    print("Saved PCA plot.")
    

def main():
    """
    Main function: cleans files, merges matrix, visualizes results.
    """
    parser = argparse.ArgumentParser(description="Process Kallisto abundance files and visualize results.")
    parser.add_argument('--input', required=True, help='Input folder containing abundance.tsv files')
    parser.add_argument('--output', required=True, help='Output folder for results')
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    abundance_files = []
    for root, dirs, files in os.walk(args.input):
        for f in files:
            if f.endswith('abundance.tsv'):
                abundance_files.append(os.path.join(root, f))
    if not abundance_files:
        print(f"[ERROR] No abundance.tsv files found in {args.input}")
        return

    cleaned_files = []
    for f in abundance_files:
        df, sample = clean_abundance_file(f, args.output)
        if df is not None:
            cleaned_files.append(os.path.join(args.output, f"{sample}_cleaned.csv"))

    merged = merge_samples(cleaned_files)
    merged.to_csv(os.path.join(args.output, 'merged_matrix.csv'))
    print(f"Merged matrix saved to {os.path.join(args.output, 'merged_matrix.csv')}")

    visualize_heatmap(merged, args.output)
    visualize_pca(merged, args.output)

if __name__ == '__main__':
    main()
