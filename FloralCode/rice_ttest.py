import pandas as pd
import numpy as np
from scipy.stats import ttest_ind

def main():
    expr_csv = 'results/rice_labeled.csv'
    out_csv = 'results/rice_ttest_results.csv'
    # Define sample-to-condition mapping
    sample_to_condition = {
        'rice_1': 'control',
        'rice_2': 'control',
        'rice_3': 'treated',
        'rice_4': 'treated'
    }
    # Load expression data
    df = pd.read_csv(expr_csv, index_col=0)
    # Get sample lists
    control_samples = [s for s, cond in sample_to_condition.items() if cond == 'control']
    treated_samples = [s for s, cond in sample_to_condition.items() if cond == 'treated']
    # Compute means
    mean_control = df[control_samples].mean(axis=1)
    mean_treated = df[treated_samples].mean(axis=1)
    # Perform t-test for each gene
    p_values = []
    for gene in df.index:
        control_vals = df.loc[gene, control_samples].values
        treated_vals = df.loc[gene, treated_samples].values
        stat, p = ttest_ind(control_vals, treated_vals, equal_var=False)
        p_values.append(p)
    # Build results DataFrame
    results = pd.DataFrame({
        'gene_id': df.index,
        'p_value': p_values,
        'mean_control': mean_control,
        'mean_treated': mean_treated
    })
    results.to_csv(out_csv, index=False)
    print(f'Results saved to {out_csv}')

if __name__ == '__main__':
    main()
