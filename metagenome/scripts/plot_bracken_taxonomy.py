import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as bpdf
import glob
import numpy as np

def get_sample_name(filepath):
    filename = os.path.basename(filepath)
    if '.pluspf' in filename:
        return filename.split('.pluspf')[0]
    return filename.split('.')[0]

def get_database_name(filepath):
    filename = os.path.basename(filepath)
    parts = filename.split('.')
    if 'bracken' in parts:
        bracken_idx = parts.index('bracken')
        if bracken_idx > 0:
            return parts[bracken_idx - 1]
    return "unknown_database"

def plot_taxonomy(suffix, output_pdf, top_n=15):
    all_files_for_suffix = glob.glob(f"results_bracken/**/*.{suffix}.tsv", recursive=True)

    if not all_files_for_suffix:
        print(f"No files found for suffix {suffix}")
        return

    database_names = sorted(list(set(get_database_name(f) for f in all_files_for_suffix)))

    for db_name in database_names:
        print(f"  Plotting for database: {db_name}, level: {suffix}...")
        
        files_for_db = [f for f in all_files_for_suffix if get_database_name(f) == db_name]

        data_frames = []
        for f in files_for_db:
            sample = get_sample_name(f)
            try:
                df = pd.read_csv(f, sep='\t')
                if 'name' in df.columns and 'fraction_total_reads' in df.columns:
                    df = df[['name', 'fraction_total_reads']]
                    df['sample'] = sample
                    data_frames.append(df)
            except Exception as e:
                print(f"Error reading {f}: {e}")

        if not data_frames:
            print(f"    No data frames generated for database {db_name}, level {suffix}. Skipping plot.")
            continue

        combined_df = pd.concat(data_frames)
        
        pivot_df = combined_df.pivot_table(index='sample', columns='name', values='fraction_total_reads', fill_value=0)
        
        if pivot_df.empty:
            print(f"    Pivot table is empty for database {db_name}, level {suffix}. Skipping plot.")
            continue

        mean_abundance = pivot_df.mean().sort_values(ascending=False)
        top_taxa = mean_abundance.head(top_n).index.tolist()
        
        plot_data = pivot_df[top_taxa].copy()
        
        remaining_cols = [col for col in pivot_df.columns if col not in top_taxa]
        if remaining_cols:
            plot_data['Other'] = pivot_df[remaining_cols].sum(axis=1)
        elif 'Other' not in plot_data.columns:
            plot_data['Other'] = 0

        plot_data = plot_data.sort_index()

        fig, ax = plt.subplots(figsize=(18, 10))
        plot_data.plot(kind='bar', stacked=True, ax=ax, width=0.8, colormap='tab20')
        
        ax.set_title(f"Taxonomic Distribution - {suffix} (Database: {db_name})", fontsize=18)
        ax.set_ylabel("Fraction of Total Reads", fontsize=14)
        ax.set_xlabel("Sample", fontsize=14)
        ax.legend(title="Taxonomy", bbox_to_anchor=(1.02, 1), loc='upper left', fontsize='small')
        plt.xticks(rotation=90, fontsize=8)
        plt.tight_layout()
        
        output_pdf.savefig(fig)
        plt.close(fig)

if __name__ == "__main__":
    pdf_filename = "taxonomic_barplots.pdf"
    with bpdf.PdfPages(pdf_filename) as pdf:
        for level in ['P', 'G', 'S']:
            print(f"Plotting level {level}...")
            plot_taxonomy(level, pdf)
    print(f"Successfully created {pdf_filename}")
