import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import os

def load_data(expr_csv, sample_to_condition):
    expr = pd.read_csv(expr_csv, index_col=0)
    cond = pd.Series(sample_to_condition)
    cond = cond.loc[expr.columns]  # Align order
    return expr, cond

def compute_log2fc(expr, cond, cond_a, cond_b):
    log_expr = np.log2(expr + 1)
    mean_a = log_expr.loc[:, cond == cond_a].mean(axis=1)
    mean_b = log_expr.loc[:, cond == cond_b].mean(axis=1)
    log2fc = mean_a - mean_b
    return log2fc

def compute_pvalues(expr, cond, cond_a, cond_b):
    log_expr = np.log2(expr + 1)
    group_a = log_expr.loc[:, cond == cond_a]
    group_b = log_expr.loc[:, cond == cond_b]
    pvals = []
    for gene in log_expr.index:
        a = group_a.loc[gene].values
        b = group_b.loc[gene].values
        if len(a) < 2 or len(b) < 2:
            pvals.append(np.nan)
        else:
            stat, p = ttest_ind(a, b, equal_var=False)
            pvals.append(p)
    return pd.Series(pvals, index=log_expr.index)

def bh_fdr(pvals):
    p = pvals.values
    n = len(p)
    idx = np.argsort(p)
    sorted_p = p[idx]
    fdr = np.empty(n)
    prev_fdr = 0
    for i in range(n-1, -1, -1):
        fdr_val = sorted_p[i] * n / (i+1)
        prev_fdr = min(fdr_val, prev_fdr) if i < n-1 else fdr_val
        fdr[i] = prev_fdr
    fdr_out = np.empty(n)
    fdr_out[idx] = fdr
    return pd.Series(fdr_out, index=pvals.index)

def plot_volcano(log2fc, pvals, fdr, cond_a, cond_b, outdir):
    plt.figure(figsize=(8,6))
    sig = fdr < 0.05
    plt.scatter(log2fc, -np.log10(pvals), c='grey', s=10, label='All genes')
    if sig.any():
        plt.scatter(log2fc[sig], -np.log10(pvals[sig]), c='red', s=10, label='Significant (FDR<0.05)')
    plt.xlabel(f'log2 Fold Change ({cond_a} vs {cond_b})')
    plt.ylabel('-log10(p-value)')
    plt.title(f'Volcano Plot: {cond_a} vs {cond_b}')
    plt.legend()
    plt.tight_layout()
    plt.show()

def main():
    expr_csv = 'results/rice_labeled.csv'  # Change as needed
    outdir = 'results'
    os.makedirs(outdir, exist_ok=True)
    sample_to_condition = {
        'rice_1': 'control',
        'rice_2': 'control',
        'rice_3': 'treated',
        'rice_4': 'treated'
        # Add more samples/conditions as needed
    }
    print("Sample-to-condition mapping:")
    for sample, cond_name in sample_to_condition.items():
        print(f"  {sample}: {cond_name}")
    expr, cond = load_data(expr_csv, sample_to_condition)
    print("Loaded conditions:", cond.unique())
    conditions = cond.unique()
    for i, cond_a in enumerate(conditions):
        for cond_b in conditions[i+1:]:
            print(f'Comparing {cond_a} vs {cond_b}...')
            log2fc = compute_log2fc(expr, cond, cond_a, cond_b)
            pvals = compute_pvalues(expr, cond, cond_a, cond_b)
            fdr = bh_fdr(pvals)
            results = pd.DataFrame({'GeneID': expr.index, 'log2FC': log2fc, 'p-value': pvals, 'FDR': fdr})
            results.to_csv(os.path.join(outdir, f'de_{cond_a}_vs_{cond_b}.csv'), index=False)
            plot_volcano(log2fc, pvals, fdr, cond_a, cond_b, outdir)
            print(f'  Results saved: de_{cond_a}_vs_{cond_b}.csv')

if __name__ == '__main__':
    main()
