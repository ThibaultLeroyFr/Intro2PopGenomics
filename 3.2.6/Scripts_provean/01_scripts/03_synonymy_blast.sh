#!/bin/bash

# Global variables
INPUT_FASTA=$1
DATABASE=$2

# Blast fasta file
cat "$INPUT_FASTA" |
    parallel -k \
    --block 1k \
    --recstart '>' \
    --pipe blastx \
    -db "$DATABASE" -query - \
    -evalue 1e-6 \
    -outfmt \"6 qseqid sseqid pident length evalue bitscore qseq sseq\" \
    -max_target_seqs 1 > \
    04_blasts/"$(basename ${INPUT_FASTA%.fasta})"."$(basename $DATABASE)"
