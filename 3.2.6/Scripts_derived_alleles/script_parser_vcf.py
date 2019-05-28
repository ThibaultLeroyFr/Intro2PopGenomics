#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
file1 = sys.argv[1] # VCF (only PASS variants)
file2 = sys.argv[2] # VCF outgroup species 1
file3 = sys.argv[3] # VCF outgroup species 2

dictsp1 = {}
dictsp2 = {}

for line2 in open(file2): # read the vcf of the first outgroup species
    if line2.startswith("#"):
        continue
    else:
        line2= line2.replace('\n','')
        splitted_line2 = line2.split('\t')
        scaffIDsp1 = splitted_line2[0] # récupérer le premier champ
        possp1 = splitted_line2[1] # récupérer le second champ
        keysp1 = (scaffIDsp1 + '-' + possp1)
        dictsp1[keysp1] = line2

for line3 in open(file3): # read the vcf of the first outgroup species
    if line3.startswith("#"):
        continue
    else:
        line3= line3.replace('\n','')
        splitted_line3 = line3.split('\t')
        scaffIDsp2 = splitted_line3[0] # récupérer le premier champ
        possp2 = splitted_line3[1] # récupérer le second champ
        keysp2 = (scaffIDsp2 + '-' + possp2)
        dictsp2[keysp2] = line3


for line in open(file1):
    line = line.replace('\n','')
    splitted_line = line.split('\t')
    scaffID = splitted_line[0] # récupérer le premier champ
    pos = splitted_line[1] # récupérer le second champ
    major = splitted_line[3] # récupérer le second champ
    minor = splitted_line[4] # récupérer le second champ
    find_pattern = scaffID + "-" + pos
    if (dictsp1.has_key(find_pattern)):
        if (dictsp2.has_key(find_pattern)):
            foundsp1 = dictsp1[find_pattern]
            foundsp2 = dictsp2[find_pattern]
            print_line = line + foundsp1 + foundsp2
            print print_line
