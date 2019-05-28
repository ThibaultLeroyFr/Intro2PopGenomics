#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
file1 = sys.argv[1] # VCF (only PASS variants)
file2 = sys.argv[2] # summaryalleles file
#1       1248    G       A       G       .       G       .
#1       1249    A       C       A       .       A       C
#1       1266    G       A       G       .       G       .
#1       1277    T       C       T       .       T       C
#1       1280    G       A       G       .       G       .
#1       1299    T       C       T       .       T       C
#1       1301    G       C       G       .       G       .
#1       1352    C       A       C       .       C       .
#1       1379    T       C       T       .       T       *

dictsp1 = {}

for line in open(file1):
    line = line.replace('\n','')
    splitted_line = line.split('\t')
    scaffID = splitted_line[0] # récupérer le premier champ
    pos = splitted_line[1] # récupérer le second
    pattern = scaffID + "-" + pos
    dictsp1[pattern] = line

# print header
print "chr\tpos\tancestral\tderived\tOb_missing_calls\tOb_dbleancestral_calls\tOb_hetero_calls\tOb_dblederived_calls\tOg_missing_calls\tOg_dbleancestral_calls\tOg_hetero_calls\tOg_dblederived_calls"

for line2 in open(file2):
    line2= line2.replace('\n','')
    splitted_line2 = line2.split('\t')
    scaffIDsp1 = splitted_line2[0] # récupérer le premier champ
    possp1 = splitted_line2[1] # récupérer le second champ
    major=splitted_line2[2]
    minor=splitted_line2[3]
    mainoutsp1 = splitted_line2[4]
    altoutsp1 = splitted_line2[5]
    mainoutsp2 = splitted_line2[6]
    altoutsp2 = splitted_line2[7]
    if altoutsp1 == altoutsp2:
        if altoutsp1 == ".": # no alternative allele for the two species
            if mainoutsp1 == mainoutsp2: # homozygous for the same allele
                # at this step we know that the two species are homozygous for the same allele
                if major == mainoutsp1: # the major allele is the ancestral
                    find_pattern = scaffIDsp1 + "-" + possp1
                    if (dictsp1.has_key(find_pattern)):
                        found_line = dictsp1[find_pattern]
                        splitted_foundline = found_line.split('\t')
                        Obarthii = splitted_foundline[9] + "\t" + splitted_foundline[10] + "\t" + splitted_foundline[11] + "\t" + splitted_foundline[12] + "\t" + splitted_foundline[13] + "\t" + splitted_foundline[14] + "\t" +  splitted_foundline[15] + "\t" + splitted_foundline[16] + "\t" + splitted_foundline[17] + "\t" + splitted_foundline[18] + "\t" +  splitted_foundline[19] + "\t" + splitted_foundline[20] + "\t" + splitted_foundline[21] + "\t" + splitted_foundline[22] + "\t" + splitted_foundline[23] + "\t" + splitted_foundline[24] + "\t" +  splitted_foundline[25] + "\t" + splitted_foundline[26] + "\t" +  splitted_foundline[52] + "\t" + splitted_foundline[53] + "\t" + splitted_foundline[54] + "\t" +  splitted_foundline[55] + "\t" + splitted_foundline[56]
                        Oglaberrima = splitted_foundline[27] + "\t" + splitted_foundline[28] + "\t" +  splitted_foundline[29] + "\t" + splitted_foundline[30] + "\t" + splitted_foundline[31] + "\t" + splitted_foundline[32] + "\t" +  splitted_foundline[33] + "\t" + splitted_foundline[34] + "\t" + splitted_foundline[35] + "\t" + splitted_foundline[36] + "\t" +  splitted_foundline[37] + "\t" + splitted_foundline[38] + "\t" + splitted_foundline[39] + "\t" + splitted_foundline[40] + "\t" + splitted_foundline[41] + "\t" + splitted_foundline[42] + "\t" +  splitted_foundline[43] + "\t" + splitted_foundline[44] + "\t" + splitted_foundline[45] + "\t" + splitted_foundline[46] + "\t" + splitted_foundline[47] + "\t" + splitted_foundline[48] + "\t" + splitted_foundline[49] + "\t" +  splitted_foundline[50] + "\t" + splitted_foundline[51]
                        missing_Obarthii = Obarthii.count("./.")
                        dble_ancestral_Obarthii = Obarthii.count("0/0")
                        heterozygous_Obarthii = Obarthii.count("0/1")
                        dble_derived_Obarthii = Obarthii.count("1/1")
                        missing_Oglaberrima = Oglaberrima.count("./.")
                        dble_ancestral_Oglaberrima = Oglaberrima.count("0/0")
                        heterozygous_Oglaberrima = Oglaberrima.count("0/1")
                        dble_derived_Oglaberrima = Oglaberrima.count("1/1")
                        to_print = scaffIDsp1 + "\t" + possp1 + "\t" + major + "\t" + minor + "\t" + str(missing_Obarthii) + "\t" + str(dble_ancestral_Obarthii) + "\t" + str(heterozygous_Obarthii) + "\t" + str(dble_derived_Obarthii) + "\t" + str(missing_Oglaberrima) + "\t" + str(dble_ancestral_Oglaberrima) + "\t" + str(heterozygous_Oglaberrima) + "\t" + str(dble_derived_Oglaberrima)
                        print to_print
                if minor == mainoutsp1: # the minor allele is the ancestral
                    find_pattern = scaffIDsp1 + "-" + possp1
                    if (dictsp1.has_key(find_pattern)):
                        found_line = dictsp1[find_pattern]
                        splitted_foundline = found_line.split('\t')
                        Obarthii = splitted_foundline[9] + "\t" + splitted_foundline[10] + "\t" + splitted_foundline[11] + "\t" + splitted_foundline[12] + "\t" + splitted_foundline[13] + "\t" + splitted_foundline[14] + "\t" +  splitted_foundline[15] + "\t" + splitted_foundline[16] + "\t" + splitted_foundline[17] + "\t" + splitted_foundline[18] + "\t" +  splitted_foundline[19] + "\t" + splitted_foundline[20] + "\t" + splitted_foundline[21] + "\t" + splitted_foundline[22] + "\t" + splitted_foundline[23] + "\t" + splitted_foundline[24] + "\t" +  splitted_foundline[25] + "\t" + splitted_foundline[26] + "\t" +  splitted_foundline[52] + "\t" + splitted_foundline[53] + "\t" + splitted_foundline[54] + "\t" +  splitted_foundline[55] + "\t" + splitted_foundline[56]
                        Oglaberrima = splitted_foundline[27] + "\t" + splitted_foundline[28] + "\t" +  splitted_foundline[29] + "\t" + splitted_foundline[30] + "\t" + splitted_foundline[31] + "\t" + splitted_foundline[32] + "\t" +  splitted_foundline[33] + "\t" + splitted_foundline[34] + "\t" + splitted_foundline[35] + "\t" + splitted_foundline[36] + "\t" +  splitted_foundline[37] + "\t" + splitted_foundline[38] + "\t" + splitted_foundline[39] + "\t" + splitted_foundline[40] + "\t" + splitted_foundline[41] + "\t" + splitted_foundline[42] + "\t" +  splitted_foundline[43] + "\t" + splitted_foundline[44] + "\t" + splitted_foundline[45] + "\t" + splitted_foundline[46] + "\t" + splitted_foundline[47] + "\t" + splitted_foundline[48] + "\t" + splitted_foundline[49] + "\t" +  splitted_foundline[50] + "\t" + splitted_foundline[51]
                        missing_Obarthii = Obarthii.count("./.")
                        dble_ancestral_Obarthii = Obarthii.count("1/1")
                        heterozygous_Obarthii = Obarthii.count("0/1")
                        dble_derived_Obarthii = Obarthii.count("0/0")
                        missing_Oglaberrima = Oglaberrima.count("./.")
                        dble_ancestral_Oglaberrima = Oglaberrima.count("1/1")
                        heterozygous_Oglaberrima = Oglaberrima.count("0/1")
                        dble_derived_Oglaberrima = Oglaberrima.count("0/0")
                        to_print = scaffIDsp1 + "\t" + possp1 + "\t" + minor + "\t" + major + "\t" + str(missing_Obarthii) + "\t" + str(dble_ancestral_Obarthii) + "\t" + str(heterozygous_Obarthii) + "\t" + str(dble_derived_Obarthii) + "\t" + str(missing_Oglaberrima) + "\t" + str(dble_ancestral_Oglaberrima) + "\t" + str(heterozygous_Oglaberrima) + "\t" + str(dble_derived_Oglaberrima)
                        print to_print

