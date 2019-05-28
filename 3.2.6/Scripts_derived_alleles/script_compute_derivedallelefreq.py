#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
file1 = sys.argv[1] # summary.counts
nb_Ob_individuals=23
nb_Og_individuals=25

#chr     pos     ancestral       derived Ob_missing_calls        Ob_dbleancestral_calls  Ob_hetero_calls        Ob_dblederived_calls    Og_missing_calls        Og_dbleancestral_calls  Og_hetero_calls        Og_dblederived_calls
#1       1248    G       A       15      0       0       8       10      0       0       15
#1       1266    G       A       15      0       0       8       10      0       0       15
#1       1280    G       A       15      0       0       8       10      0       0       15
#1       1301    G       C       15      0       0       8       10      0       0       15
#1       1352    C       A       0       21      0       2       0       25      0       0

for line in open(file1):
    line = line.replace('\n','')
    splitted_line = line.split('\t')
    chrom=splitted_line[0]
    if chrom != "chr": # all lines except the header
        pos=splitted_line[1]
        ancestral=splitted_line[2]
        derived=splitted_line[3]
        missingOb=int(splitted_line[4])
        missingOg=int(splitted_line[8])
        missingObrate=float(missingOb)/float(nb_Ob_individuals)
        missingOgrate=float(missingOg)/float(nb_Og_individuals)
        if missingObrate < 0.1 and missingOgrate < 0.1: # no more than 10% of missing data for each individual
            derivedallelefreqOb=((float(splitted_line[7])*2)+float(splitted_line[6]))/(float(nb_Ob_individuals)*2)
            derivedallelefreqOg=((float(splitted_line[11])*2)+float(splitted_line[10]))/(float(nb_Og_individuals)*2)
            print chrom + "\t" + pos + "\t" + ancestral + "\t" + derived + "\t" + str(missingObrate) + "\t" + str(missingOgrate) + "\t" + str(derivedallelefreqOb) + "\t" + str(derivedallelefreqOg)
