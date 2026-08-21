import os
import glob
import json

# Initialize dictionary of gene families
db = {'PGM/TPB2':[], 'TPB1':[], 'TPB6':[]}

# Check Eukphylo database for fasta files with the orthogroups in the dictionary and add respective taxa to the dictionary
for f in glob.glob('/sciclone/scr10/abharathan01/transposases/DB_FASTA/Eukphylo_AAs/*fasta'):
	taxon=f.rpartition("/")[-1].partition(".")[0]
	if open(f).read().count('OG6_101171') > 0:
                db['PGM/TPB2'].append(taxon)
	if open(f).read().count('OG6_102536') > 0:
                db['TPB1'].append(taxon)
	if open(f).read().count('OG6_111972') > 0 or open(f).read().count('OG6_100364') > 0 or open(f).read().count('OG6_102536') > 0:
                db['TPB6'].append(taxon)


	
# Make upset plot using updated dictionary
from upsetplot import UpSet
from upsetplot import from_contents
import matplotlib.pyplot as plt
test_data = from_contents(db)

#Make taxon list
with open("PiggyBacs.json", "w", encoding = "utf-8") as f:
        json.dump(db, f, indent=4)

upset_plot = UpSet(test_data, subset_size = "count", sort_by='cardinality').plot()
plt.title("Ciliate PiggyBacs")
plt.savefig("Ciliate_PiggyBacs_UpsetPlot.png", dpi=300)
