#!/bin/bash

input_file=$1
input_basename=$(basename "$input_file")
path="/home/qurou/software/provean-1.1.5/"

output_scores="08_provean_scores"
output_sseq="07_provean_sseq"
output_sss="$output_sseq"/${input_basename%.fasta}.sss

"$path"/provean.sh \
    -q $input_file \
    -v ${input_file%.fasta}.var \
    --save_supporting_set "$output_sss"  >> "$output_score"/allout${input_basename%.fasta}.txt 2>&1 
