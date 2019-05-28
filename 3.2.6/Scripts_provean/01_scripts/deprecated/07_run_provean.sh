#!/bin/bash

#ncpus=$1

TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

#if [[ -z "$ncpus" ]];
#then
#    ncpus=10
#fi

path="/path/to/software/provean-1.1.5/"

input="06_provean"
output="07_provean_scores"

if [[ ! -d "$output" ]]
then
    echo "creating output dir"
    mkdir "$output"
fi

for fast in $(ls -1 $input/*fasta )
do
    echo "analysing $fast"

    name=$(basename $fast )

    $path/provean.sh \
    -q $input/$name \
    -v $input/${name%.fasta}.var  --save_supporting_set $output/${name%.fasta}.sss \
    '2>&1' '>>' 10-log_files/"$TIMESTAMP"_provean"${name%.fasta}".log
done 
