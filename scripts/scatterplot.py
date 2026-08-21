import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os 
import glob

#x axis is percent identity, y axis is bit score

folder_path = "../ciliate_TSVs/Controls_PiggyBacs/*.tsv"

df_list = []
for f in glob.glob(folder_path):
        if os.path.getsize(f) == 0:
                print(f"empty file: {os.path.basename(f)}")
                continue
        df = pd.read_csv(f, sep="\t", header=None)
        clean = os.path.basename(f)
        clean = clean.replace("Protein-Matching-", "")
        clean = clean.replace(".Test.EukPhylo.tsv", "")
        length = df[9] - df[8]
        df["source_file"] = clean
        df["Percent_Identity"] = df[2]
        df["Bit_Score"] = df[11]/length
#        df_filtered = df[(df['Percent_Identity'] >= 20) & (df['Percent_Identity'] <= 100)]
        df_list.append(df)

                 
 
combined_df = pd.concat(df_list, ignore_index=True)


scatplot = sns.scatterplot(data=combined_df, x="Percent_Identity", y="Bit_Score", hue="source_file", alpha=1)


plt.legend(bbox_to_anchor = (1.05, 1), loc="upper left", fontsize = "x-small")
plt.xlabel("Percent Identity")
plt.ylabel("Normalized Bit Score")
plt.title("Eukaryotic PiggyBac Alignment Performance")
plt.tight_layout()
plt.savefig("PiggyBac_v_Control_Euk_Normalized_Bit_Score_Scatterplot.png", dpi=300)
scatplot.set_xlim(0,100)
scatplot.set_xticks(np.arange(0,101,20))
