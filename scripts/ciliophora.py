import os
import csv
from collections import Counter
file_path = "../DB_FASTA/Phylogenomic_Datasets.csv"
phlist = []

#Count ciliate taxa groups
if os.path.exists(file_path):
    with open(file_path, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        next(reader, None)
        for row in reader:
            if len(row) > 8:
                phylum = row[8]
                if phylum == "Ciliophora":
                    phlist.append(phylum)

ciliate = Counter(phlist)
print(ciliate)
