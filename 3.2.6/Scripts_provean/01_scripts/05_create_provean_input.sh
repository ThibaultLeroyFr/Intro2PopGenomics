#!/bin/bash

input=$1

if [[ -z "$input" ]]
then 
    echo "error you must provide input name "
    echo "input is obtained from 01_scripts/04_find_synonymy_with_filters.py"
fi

awk '{print $1"\t"$4"\t"$7 }' $input | sed '1d' > ${input}.for_provean.tsv
