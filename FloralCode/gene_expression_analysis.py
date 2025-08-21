import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
import os

# --- Load expression matrix and condition labels ---
def load_data(expr_csv, sample_to_condition):
    expr = pd.read_csv(expr_csv, index_col=0)
    cond = pd.Series(sample_to_condition)
    cond = cond.loc[expr.columns]  # align conditions to samples
    return expr, cond

# --- Compute log2 fold change between two conditions ---
def compute_log2fc(expr, cond, cond_a, cond_b):
    log_expr = np.log2(expr + 1)
    mean_a = log_expr.loc[:, cond == cond_a].mean(axis=1)
    mean_b = log_expr.loc[:, cond == cond_b].mean(axis=1)
    return mean_a - mean_b

# --- Compute p-values (simple t-test here, easier than permutation for all pairs) ---
def compute_pvalues(expr, cond, cond_a, cond_b):
    log_expr = np.log2(expr + 1)
    pvals = []
    for gene in log_expr.index:
        vals_a = log_expr.loc[gene, cond == cond_a].values
        vals_b = log_expr.loc[gene, cond == cond_b].values
        if len(vals_a) < 2 or len(vals_b) < 2:
            pvals.append(np.nan)
        else:
            stat, p = ttest_ind(vals_a, vals_b, equal_var=False)
            pvals.append(p)
    return pd.Series(pvals, index=log_expr.index)

# --- FDR correction (Benjamini–Hochberg) ---
def bh_fdr(pvals):
    p = pvals.fillna(1).values
    n = len(p)
    idx = np.argsort(p)
    sorted_p = p[idx]
    fdr = np.empty(n)
    prev_fdr = 1.0
    for i in range(n-1, -1, -1):
        fdr_val = sorted_p[i] * n / (i+1)
        prev_fdr = min(fdr_val, prev_fdr)
        fdr[i] = prev_fdr
    fdr_out = np.empty(n)
    fdr_out[idx] = fdr
    return pd.Series(fdr_out, index=pvals.index)

# --- Volcano plot ---
def plot_volcano(log2fc, pvals, fdr, cond_a, cond_b, outdir, fdr_thresh=0.05):
    plt.figure(figsize=(8,6))
    plt.scatter(log2fc, -np.log10(pvals), c='grey', s=10, label='All genes')
    sig = fdr < fdr_thresh
    if sig.any():
        plt.scatter(log2fc[sig], -np.log10(pvals[sig]), c='red', s=10, label=f'Significant (FDR<{fdr_thresh})')
    plt.xlabel(f'log2 Fold Change ({cond_a} vs {cond_b})')
    plt.ylabel('-log10(p-value)')
    plt.title(f'Volcano Plot: {cond_a} vs {cond_b}')
    plt.legend()
    os.makedirs(outdir, exist_ok=True)
    plt.savefig(os.path.join(outdir, f'volcano_{cond_a}_vs_{cond_b}.png'), bbox_inches='tight')
    plt.close()

# --- Main workflow ---
def main():
    expr_csv = 'results/rice_labeled.csv'
    expr = pd.read_csv(expr_csv, index_col=0)
    sample_names = list(expr.columns)

    # Example mapping: assign each sample to one of many conditions
    all_conditions = ['Heat','Water','Osmosis','Cond4','Cond5','Cond6','Cond7','Cond8']
    sample_to_condition = {sn: all_conditions[i % len(all_conditions)] for i, sn in enumerate(sample_names)}

    expr, cond = load_data(expr_csv, sample_to_condition)

    # Loop over all unique condition pairs
    conditions = cond.unique()
    outdir = "results/pairwise"
    os.makedirs(outdir, exist_ok=True)

    for i in range(len(conditions)):
        for j in range(i+1, len(conditions)):
            cond_a, cond_b = conditions[i], conditions[j]
            print(f"Comparing {cond_a} vs {cond_b}...")

            log2fc = compute_log2fc(expr, cond, cond_a, cond_b)
            pvals = compute_pvalues(expr, cond, cond_a, cond_b)
            fdr = bh_fdr(pvals)

            results = pd.DataFrame({
                'GeneID': expr.index,
                f'log2FC_{cond_a}_vs_{cond_b}': log2fc,
                'p-value': pvals,
                'FDR': fdr
            })

            results.to_csv(os.path.join(outdir, f'{cond_a}_vs_{cond_b}.csv'), index=False)

            plot_volcano(log2fc, pvals, fdr, cond_a, cond_b, outdir)

if __name__ == '__main__':
    main()
