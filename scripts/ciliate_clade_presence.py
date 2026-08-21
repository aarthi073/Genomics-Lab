import csv
import os
import sys
from collections import defaultdict
from collections import Counter
import matplotlib.pyplot as plt
import pandas as pd
import json

file_path = "../DB_FASTA/Phylogenomic_Datasets.csv"

#pass json file of taxa that have the gene family
try:
    j_path = sys.argv[1]
except IndexError:
    print("no database passed")
    sys.exit()

with open(j_path) as f:
    db = json.load(f)

#replace "_" with a space in the taxa list
clean_to_orig = {}
for species_list in db.values():
    for species in species_list:
        clean_name = species.replace('_', ' ')
        clean_to_orig[clean_name] = species

clade_results = defaultdict(lambda: defaultdict(list))
phlist_global = [] # every (major, minor) pair from the phylogenomic dataset
total_global = 0 #total number of taxa in phylogenomic dataset
seen_pairs = set() #Unique (major, minor) pairs

if os.path.exists(file_path):
    with open(file_path, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        next(reader, None)
        for row in reader:
            total_global += 1
            taxon = row[0]  
            major_clade = row[7]
            minor_clade = row[8]
            pair = (major_clade, minor_clade)
            
            if pair not in seen_pairs:
                seen_pairs.add(pair)
               # phlist_global.append(pair)
           
            #match taxa and (major clade, minor clade) between phylogenomic dataset and taxa list for gene family
            if taxon in clean_to_orig:
                orig_species_name = clean_to_orig[taxon]
                clade_pair = [major_clade, minor_clade]
                                
            #append clade pair to clade_results if taxa match between the phylogenomic dataset and taxa list for each respective gene family
                for gene_family, species_list in db.items():
                    if orig_species_name in species_list:
                        if clade_pair not in clade_results[gene_family][orig_species_name]:
                            clade_results[gene_family][orig_species_name].append(clade_pair)


# Write one output file per gene family
total_species_all_families = sum(len(species_list) for species_list in db.values())
for gene_family in clade_results:
    output_tsv = f"../Docs/{gene_family}_clade_results.tsv"
    pct_cov = len(clade_results[gene_family])/total_species_all_families
    # Collect all clade pairs for THIS gene family only
    pairs_this_family = []
    for species_name, clade_pair in clade_results[gene_family].items():
        for major_clade, minor_clade in clade_pair:
            pairs_this_family.append((major_clade, minor_clade))
            
    total_this_family = len(pairs_this_family)
    pair_counts_this_family = Counter(pairs_this_family)
    Percentage_this_family = {
        pair: (count / total_this_family) * 100 for pair, count in pair_counts_this_family.items()
    }
    # Build table rows with both percentages
    table_data = []
    for pair_key, pct_family in Percentage_this_family.items():
        major_clade, minor_clade = pair_key
        
        table_data.append({
            "Major_Clade": major_clade,
            "Minor_Clade": minor_clade,
            "Percentage_in_Taxonomic_Group": round(pct_family, 2),
            "Percentage_in_Dataset": round(pct_cov, 2)
        })

   # Built tables        
    with open(output_tsv, "w", newline="", encoding="utf-8") as f:
        fieldnames = ["Major_Clade", "Minor_Clade", "Percentage_in_Taxonomic_Group", "Percentage_in_Dataset"]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(table_data)
        
    # Stacked barplot
    df = pd.read_csv(output_tsv, sep="\t")
    df = df.set_index(["Major_Clade","Minor_Clade"])
    df = df[["Percentage_in_Taxonomic_Group"]].fillna(0)
    df.plot(kind="barh", figsize=(25,10))
    plt.title("Gene Family Presence in Each Taxonomic Group", fontsize=20)
    plt.xlabel("Percentage in Taxonomic Group", fontsize=20)
    plt.ylabel("Taxonomic Group", fontsize=20)   
    plt.savefig("../Figures/Gene_Family_Presence_Barplot.png", dpi=300)
