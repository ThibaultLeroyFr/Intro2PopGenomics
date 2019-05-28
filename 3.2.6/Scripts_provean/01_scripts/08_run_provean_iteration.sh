#!/bin/bash

input_file=$1
chr=$2

input_basename=$(basename "$input_file")
#path="/home/qurou/software/provean-1.1.5/"
#path="/home/quentin/scratch/24.rice/provean-1.1.5/"
path="/home/quentin/software/provean-1.1.5/"

output_folder="07_provean_match"
output_sss="$output_folder"/"$chr"/${input_basename%.fasta}.sss

"$path"/provean.sh \
    -q $input_file \
    -v ${input_file%.fasta}.var \
    --save_supporting_set "$output_sss"  >> 08_provean_scores/"$chr"/score_${input_basename%.fasta}.txt 2>&1 
