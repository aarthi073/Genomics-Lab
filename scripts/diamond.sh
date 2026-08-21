#!/bin/bash
#SBATCH --job-name=diamond_test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH -t 00:60:00
#SBATCH --mem=12G



#compares protein sequences against reference protein database (Eukphylo); alternative to NCBI blastp
diamond blastp \
	-k 1 \
	-f 6 \
	-o ../ciliate_TSVs/B.stoltei_5.EukPhylo.tsv \
	-q ../query_fasta/original_ciliate/B.stoltei_5.fasta\
        -d ../DB_FASTA/eukphylo_db.updated.no_inframe_stop.fasta\
	--ultra-sensitive 
