#Allows for interaction with operating system 
import os
import csv
from collections import Counter

#need to keep updating this depending on what TSV directory you are working with.
folder_path = "../ciliate_TSVs"

#Initialize empty list for table data 
table_data = []

#For each file in the current directory...
for filename in os.listdir(folder_path):

#That ends with tsv... 
    if filename.endswith(".tsv"):
	
	#Initialize empty list to store all orthogroups for each gene 
        OG6 = []
	
#creates a file path that combines the current folder name with the file name to use the correct
#path naming syntax depending on whatever operating system you are working on.
        file_path = os.path.join(folder_path, filename)

        with open(file_path, "r") as f:

#read the tsv file
            reader = csv.reader(f, delimiter="\t")


            #next(reader, None)  


            for row in reader:
                try:
                    e_value = float(row[10])

#get the OG number after |

                    gene_family = row[1].split("|")[-1]


                    # filter hits based on e-value

                    if e_value < 1e-20 and "OG" in gene_family:
                        OG6.append(gene_family)

                except (ValueError, IndexError):
                    continue

        gene_counts = Counter(OG6)

#OG with most frequency

        if gene_counts:
            #top_og, count = gene_counts.most_common(3)[0]
            top_og, count = gene_counts.most_common(1)[0]
        #sorted_counts = list(gene_counts.items())
        #sorted_counts.sort(key=lambda x:x[1],reverse=True)

        else:
            top_og, count = "None", 0
        
        sorted_counts = list(gene_counts.items())
        sorted_counts.sort(key=lambda x:x[1],reverse=True)



        if sorted_counts:
            print(filename, sorted_counts[:3])
        else:
            print(filename, "no orthogroups found")


        
        top1 = sorted_counts[0] if len(sorted_counts) > 0 else ("no orthogroup", 0)
        top2 = sorted_counts[1] if len(sorted_counts) > 1 else ("no orthogroup", 0)
        top3 = sorted_counts[2] if len(sorted_counts) > 2 else ("no orthogroup", 0)
        top4 = sorted_counts[3] if len(sorted_counts) > 3 else ("no orthogroup", 0)

        table_data.append({
            "Filename": filename,
            "Top_Gene_Family": top1[0],
            "Top_Count": top1[1], 
            "2_Gene_Family": top2[0],
            "2_Count": top2[1],
            "3_Gene_Family": top3[0],
            "3_Count": top3[1],
            "4_Gene_Family": top4[0],
            "4_Count": top4[1]
        })


output_file = "../ciliate_TSVs/summary_table.tsv"

#Write summary table
with open(output_file, "w", newline="") as f:
    fieldnames = ["Filename", "Top_Gene_Family", "Top_Count", "2_Gene_Family", "2_Count", "3_Gene_Family", "3_Count", "4_Gene_Family", "4_Count"]
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")

    writer.writeheader()
    writer.writerows(table_data)

print(f"Saved: {output_file}")
