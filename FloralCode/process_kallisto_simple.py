import os
import pandas as pd
import numpy as np

# ---------------------------
# LOAD ALL THE FILES AND HOPE THEY EXIST
# ---------------------------
# This function tries to find every 'abundance*.tsv' file in the folder for a species.
# Each file is a sample. Each row is a gene. Every file screams at us with numbers.
# It reads them, renames the count column with the sample name, sets gene IDs as index,
# and merges everything into one messy-but-usable matrix.
# Returns both the expression matrix and a list of sample names because apparently we need both.
def load_expression_matrix(species_folder):
    dfs = []
    col_names = []
    found_files = []
    for root, dirs, files in os.walk(species_folder):
        for f in files:
            if f.startswith('abundance') and f.endswith('.tsv'):
                path = os.path.join(root, f)
                found_files.append(path)
                col_name = os.path.splitext(f)[0]
                df = pd.read_csv(path, sep='\t', usecols=['target_id', 'est_counts'])
                df = df.rename(columns={'est_counts': col_name}).set_index('target_id')
                dfs.append(df)
                col_names.append(col_name)
    print(f"  Found {len(found_files)} abundance*.tsv files in {species_folder}:")
    for fp in found_files:
        print(f"    {fp}")
    if not dfs:
        print(f"  [ERROR] No abundance*.tsv files found in {species_folder}")
        return None, []
    # Merge all the samples together and fill missing data with zeros
    expr = pd.concat(dfs, axis=1).fillna(0)
    expr = expr[sorted(expr.columns)]  # Sort columns just to make things orderly
    return expr, col_names

# ---------------------------
# TURN NUMBERS INTO SOMETHING WE CAN ACTUALLY PLOT
# ---------------------------
# Log2 transform to calm down wildly expressed genes
# Z-score normalize so every gene behaves nicely
# Basically makes the data look like we know what we’re doing
def preprocess_matrix(expr):
    expr = np.log2(expr + 1)
    expr = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    return expr

# ---------------------------
# MAIN FUNCTION: TRY NOT TO DIE
# ---------------------------
# Loops through all species we care about (rose and rice for now)
# Loads their expression matrices, cleans them, and saves them
# Computes correlation matrices for samples, because someone said "coexpression"
# Saves the results in 'output' folder.
# Prints warnings if files are missing or if there aren’t enough samples
# Finally writes a summary.txt so we can pretend we’re organized
def main():
    species_folders = {
        'rose': 'data/rose',
        'rice': 'data/rice',
        # Add more species as needed
    }
    outdir = 'output'
    os.makedirs(outdir, exist_ok=True)
    summary_lines = []
    for species, folder in species_folders.items():
        print(f"Processing {species}...")
        expr, col_names = load_expression_matrix(folder)
        if expr is None or expr.shape[1] == 0:
            summary_lines.append(f"{species}: 0 genes/transcripts, 0 samples (no files found)")
            print(f"[WARN] Skipping {species} due to no data.")
            continue
        summary_lines.append(f"{species}: {expr.shape[0]} genes/transcripts, {len(col_names)} samples")
        expr = preprocess_matrix(expr)
        if expr.shape[1] < 2:
            print(f"[WARN] Not enough samples to compute correlation for {species} (found {expr.shape[1]} sample(s))")
            continue
        # Compute Pearson correlation between samples because why not
        corr = expr.T.corr(method='pearson')
        corr.to_csv(os.path.join(outdir, f'{species}_correlation.csv'))
        print(f"Saved {species}_correlation.csv")
        print("Sample correlation values (first 5x5):")
        print(corr.iloc[:5, :5])
    with open(os.path.join(outdir, 'summary.txt'), 'w') as f:
        f.write('\n'.join(summary_lines))
    print("Summary written to summary.txt")

# ---------------------------
# RUN EVERYTHING AND PRAY IT WORKS
# ---------------------------
if __name__ == '__main__':
    main()
