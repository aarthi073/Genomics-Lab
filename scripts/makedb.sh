#!/bin/bash
#SBATCH --job-name=blast_makedb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH -t 00:60:00
#SBATCH --mem=12G

#blast germline database to output database formatted for blastn
DB="/sciclone/scr10/abharathan01/transposases/DB_FASTA"

makeblastdb -in ${DB}/Tetrahymena_thermophila_mic_assembly_v6.fasta -dbtype nucl -parse_seqids -out T.therm_mic_db
