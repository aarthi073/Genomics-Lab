import json
import os

#Define file paths
file_path = "./PiggyBac_PBLE_Taxa.json"

output = os.path.expanduser("./PiggyBac_PBLE_Trees.taxon_list.txt")

#load the JSON data and add taxa to a set to reject duplicates
taxa = set()
with open(file_path, "r") as f:
	data = json.load(f)
	for i in data.values():
		for v in i:
			taxa.add(v)

# open the output file and write items

with open(output, "w") as of:
	#grab values
	for v in sorted(taxa):
		of.write(f"{v}\n")



