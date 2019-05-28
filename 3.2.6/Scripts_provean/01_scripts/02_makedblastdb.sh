#!/bin/bash

input=$1 #protein fasta file stored in 02_data/ folder
title=$2 #title for the blast db

if [[ -z $input ]]
then
    echo "error input protein file needed!!!"
fi

if [[ -z "$title" ]]
then
    echo "error alias for the database needed"
fi    

input_type="fasta"
dbtype="prot"

makeblastdb -in "$input" \
            -input_type "$input_type" \
            -dbtype "$dbtype"\
            -title "$title"
