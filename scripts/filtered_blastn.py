from collections import defaultdict
import sys
import itertools
import csv
import pandas as pd
import glob
import os




def main():
        try:
                domtbl_files = sys.argv[1:]
        except IndexError:
                print("no domain table passed")
                sys.exit()
        first_dict=defaultdict(list)
        second_dict={}
        third_dict={}
        fourth_dict = {}

        headers = ["qseqid", "sseqid", "pid","length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

        #convert txt file to dataframe
        def dataframe(domtbl_file):
                cols = ["qseqid", "sseqid", "pid","length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
                rows = []
                with open(domtbl_file) as f:
                        for line in f:
                                if line.startswith("#"):
                                    continue
                #remove trailing white space and treat whitespace as a single separator. without empty entries where there are
                #more than one space, use None to avoid empty entries. 11 is the max number of splits.
                                fields = line.strip().split(None, 11)
                                rows.append(fields)
                df = pd.DataFrame(rows, columns=cols)

                numeric_cols = ["pid","length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"] 
                   
                for col in numeric_cols:
                    df[col] = pd.to_numeric(df[col], errors="coerce")
                df.to_csv(f"../MIC/blastn/df/{name}_df.tsv", sep="\t", index=False)
                return df
        
        def e_value(domtbl):
                domtbl.columns=domtbl.columns.str.strip()

                evalue = domtbl[domtbl["evalue"] < 1e-15]
                for seq, val in evalue.groupby("qseqid"):
                        first_dict[seq] = [tuple(row) for row in val.values]
                return first_dict

        def norm(d1):
        #Low
                for k, v in d1.items():
                        i_list = []
                        for i in v:
                                subj_length = int(i[9]) - int(i[8])
                                if subj_length != 0:
                                        normalized_score = round(i[11] / subj_length, 2)
                                else:
                                        normalized_score = 0.0
                                i_list.append(i + (normalized_score,))
                        if i_list:
                                second_dict[k] = i_list
                        for k in second_dict:
                                second_dict[k] = sorted(second_dict[k], key=lambda i: i[-1], reverse=True)
                return second_dict

        # extract entries with a bit score that is greater than 1 
        def bit(d2):
                # Medium
                for k,v in d2.items():
                        i_list = []
                        for i in v:
                                if i[-1] > 1.0:
                                        i_list.append(i)
                # append i tuple that has coverage greater than 0.4 for that key
                        if i_list:
                                third_dict[k] = i_list

                return third_dict




       # Extract entries with alignments that only overlap less than 40%
        def overlap(d3):
                # High
                for k,v in d3.items():
                        kept = []
                        #Takes care of if there is only one tuple for a sequence
                        if len(v) == 1:
                                for i in v:
                                        kept.append(i)
                        best_score =  max(itertools.chain.from_iterable(d3.values()), key=lambda i: i[-1])
                        best_from, best_to = int(best_score[8]), int(best_score[9])
                        best_len = best_to-best_from
                        for i in v:
                                if i == best_score:
                                        kept.append(i)
                                        continue
                                if i[9] - i[8] < 50:
                                        continue
                                i_from, i_to = int(i[8]), int(i[9])
                                overlap_start = max(i_from, best_from)
                                overlap_end = min(i_to, best_to)
                                overlap_len = max(0, overlap_end - overlap_start)
                                if overlap_len <= (0.4 * best_len):
                                        kept.append(i)
                        if kept:
                                fourth_dict[k] = kept
                        
                return fourth_dict
        #only called when open reading frame translation results are analyzed        
        def best_domain(d4):
                best = max(itertools.chain.from_iterable(d4.values()), key=lambda i: i[-1])
                print(f"best across all reading frames: {best}") 

                return best

        def process_reading_frame(domtbl_subset, rf, sample):
                if domtbl_subset.empty:
                        return None
               
                d1 = e_value(domtbl_subset)
                d2 = norm(d1)
                d3 = coverage(d2)
                d4 = overlap(d3)
                d5 = best_domain(d4)
                return d5

        def rf_best(domtbl_files):
                all_results = {}
                for domtbl_file in domtbl_files:
                    filename = os.path.basename(domtbl_file)
                    sample_name = "_".join(filename.split("_")[:2])
                    
                    #load dataframe
                    domtbl = dataframe(domtbl_file)
                    
                    #data for each reading frame
                    Rf1_data = domtbl[domtbl['Sequence'].str.contains('Rf1', na=False)].copy()
                    Rf2_data = domtbl[domtbl['Sequence'].str.contains('Rf2', na=False)].copy()
                    Rf3_data = domtbl[domtbl['Sequence'].str.contains('Rf3', na=False)].copy()
                    Rf4_data = domtbl[domtbl['Sequence'].str.contains('Rf4', na=False)].copy()
                    Rf5_data = domtbl[domtbl['Sequence'].str.contains('Rf5', na=False)].copy()
                    Rf6_data = domtbl[domtbl['Sequence'].str.contains('Rf6', na=False)].copy()
   
                 
                    Rf1_results = process_reading_frame(Rf1_data, "Rf1", sample_name)
                    Rf2_results = process_reading_frame(Rf2_data, "Rf2", sample_name)
                    Rf3_results = process_reading_frame(Rf3_data, "Rf3", sample_name)
                    Rf4_results = process_reading_frame(Rf4_data, "Rf4", sample_name)
                    Rf5_results = process_reading_frame(Rf5_data, "Rf5", sample_name)
                    Rf6_results = process_reading_frame(Rf6_data, "Rf6", sample_name)
                    
                    all_results[sample_name] = {
                        'Rf1': Rf1_results,
                        'Rf2': Rf2_results,
                        'Rf3': Rf3_results,
                        'Rf4': Rf4_results,
                        'Rf5': Rf5_results,
                        'Rf6': Rf6_results
                    }
                
                return all_results



                                
        #keep track of which dictionary you are passing into the table when you make an argume>
        def make_table(d4):
                rows = []
                for k,v in d4.items():
                        rows += list(set(v))

                df = pd.DataFrame(rows, columns = headers + ["Normalized_Score"])
                df.to_csv(f"../MIC/blastn/Bit/{name}_Bit_Table.tsv", sep="\t", index=False)
                return df


        
        for domtbl_file in domtbl_files:
                filename = os.path.basename(domtbl_file)
                name = "_".join(filename.split("_")[:2])
                domtbl = dataframe(domtbl_file)
                d1 = e_value(domtbl)
                d2 = norm(d1)
                d3 = bit(d2)
                d5 = make_table(d3)

#        print(rf_best(domtbl_files)) 
       # print(domain_per_rf)               
main()
