#!/bin/bash
#SBATCH --job-name=tblastn
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH -t 00:60:00
#SBATCH --mem=12G


#Compare germline PiggyBac nucleotide query sequences to find similar regions.

DB="/sciclone/scr10/abharathan01/transposases/DB_FASTA"

blastn -query T.therm_PiggyBac_library.fasta -db ${DB}/T.therm_mic_db -evalue 1e-20 -out blastn_results.txt -outfmt 6
