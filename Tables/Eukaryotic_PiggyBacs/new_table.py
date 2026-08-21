import csv
import os
import glob
from collections import defaultdict
from collections import Counter
import pandas as pd


folder_path = "/sciclone/scr10/abharathan01/transposases/PiggyBac_TSVs/*_results.tsv"

family = []


db  = {'OG6_103189': ["Euastrum_humerosum",
        "Serratifera_varisterna",
        "Nitzschia_laevis"], 'OG6_105171': ["Montipora_capricornis",
        "Zygnemopsis_sp",
        "Gallus_gallus"], 'OG6_101171': ["Monocystis_agilis",
        "Polypterus_senegalus",
        "Lithophyllum_stictiforme_cabiochiae"]}


for file_path in glob.glob(folder_path):
    file_name = os.path.basename(file_path)
    gene_family_name = os.path.splitext(file_name)[0]

    df = pd.read_csv(file_path, sep="\t")
    df["Gene_Family"] = gene_family_name
    family.append(df)

    master_df = pd.concat(family, ignore_index = True)

    taxon_counts = (
        master_df.groupby(["Major_Clade", "Taxon"])["Gene_Family"]
        .nunique()
        .reset_index()
    )
    taxon_counts.columns = ["Major_Clade", "Taxon", "GF_Count"]

    summary_table = (
    taxon_counts.groupby("Major_Clade")["GF_Count"]
    .agg(
        Unique_Taxa="count",  # Total unique species in this clade
        Avg_GF_per_Taxon="mean",  # Average number of families those taxa are in
        Taxa_in_More_Than_One_GF=lambda x: (x > 1).sum(),  # Count of taxa in >1 family
    )
    .reset_index()
    )

    summary_table.to_csv("major_clades_summary.tsv", sep="\t", index=False)
    print(summary_table)


