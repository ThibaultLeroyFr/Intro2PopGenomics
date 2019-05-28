#!/bin/bash

#script to explore length and seq identity 
#to choose a set of parameters 
#to feed 01_scripts/03_find_synonymy_with_filters.py 

input=$1
if [[ -z "input" ]]
then
    echo "error must provide input name!!! "
fi

awk '{print $3, $4}' "$input" > 04_blasts/data_to_plot

Rscript 01_scripts/utility_scripts/ploting_data.R
