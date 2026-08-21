#!/bin/bash
#SBATCH --job-name=blastx
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH -t 00:60:00
#SBATCH --mem=12G



#compares coding sequences against reference tetrahymena protein database (Eukphylo); alternative to NCBI blastx
diamond blastx \
        -k 1 \
        -f 6 \
        -o tetra_pgbd_TSVs/Protein-Matching-Tetrahymena_vorax.Test.EukPhylo.tsv \
        -q query_fasta/original_tetra_pgbd/Tetrahymena_vorax_coding_sequences.fasta.fasta\
        -d DB_FASTA/Tetrahymena_pgbd.fasta\
        --ultra-sensitive
