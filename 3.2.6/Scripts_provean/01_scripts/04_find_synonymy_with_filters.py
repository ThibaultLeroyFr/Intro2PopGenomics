#!/usr/bin/env python
"""From blast result in special format (see below) find hits for both variants
and pair them to find synonymy.

Usage:
    ./01_scripts/04_find_synonymy_with_filters.py input_file min_simil min_len output_file

Blast format:
    100054_28_C 01350001 80.000 25 3.36e-07 42.7 QRVVYMAQYITGTKLPAIQELYTRQ QRVVRTAQYITGAKLPAIQDLYNRQ
    100054_28_C 73373001 80.000 25 1.65e-06 42.0 QRVVYMAQYITGTKLPAIQELYTRQ QRVVRTAQYITGAKLPAIQDLYTRR
    100054_28_C 06681001 76.000 25 2.44e-06 40.4 QRVVYMAQYITGTKLPAIQELYTRQ QRVVRTAQYITGAKLPAIQDFYTRR
    100054_28_C 78925001 76.000 25 2.47e-06 40.0 QRVVYMAQYITGTKLPAIQELYTRQ QRVVRTAQYITGATLPAIQNLYTRR
"""

# Modules
from collections import defaultdict
import sys

# Parse user options
try:
    input_file = sys.argv[1]
    min_simil = float(sys.argv[2])
    min_len = int(sys.argv[3])
    output_file = sys.argv[4]
except:
    print __doc__
    sys.exit(1)

# Put input_file in dictionary
locus_dict = defaultdict(lambda: defaultdict(list))

with open(input_file) as infile:
    for line in infile:
        query, hit, simil, length, evalue, score, qseq, sseq  = line.strip().split()
        length = int(length)
        simil = float(simil)
        info = [query, hit, simil, length, evalue, score, qseq, sseq]

        if (length < min_len):
            continue

        if (simil < min_simil):
            continue

        locus = "_".join(query.split("_")[0:3])
        locus_dict[locus][query].append(info)

# Get best matches
with open(output_file, "w") as outf:
    outf.write("LocusID1\tPosition\tProteinID\tAlignment1\tevalue1\tLocusID2\tAlignment2\tevalue2\tSameLength\tSynonymy\n")
    for locus in locus_dict:
        best_hits = []
        best_evalue = 1.0

        if len(locus_dict[locus]) == 2:
            #print "good_locus"
            v1, v2 = locus_dict[locus].keys()

            # Create pairs of hits to the same protein
            for hit1 in locus_dict[locus][v1]:
                gene1 = hit1[1]
                evalue1 = float(hit1[4])

                for hit2 in locus_dict[locus][v2]:
                    gene2 = hit2[1]

                    if gene2 == gene1:
                        evalue2 = float(hit2[4])
                        lowest_evalue = min([evalue1, evalue2])

                        if lowest_evalue < best_evalue:
                            best_evalue = lowest_evalue
                            best_hits = [hit1, hit2]

            #print best_hits
            if best_hits:
                b1, b2 = best_hits

                same_length = 0
                if len(b1[6]) == len(b2[6]):
                    same_length = 1

                synonymy = "1"
                if b1[6] != b2[6]:
                    synonymy = "0"

                position = b1[0].split("_")[2]
                line = b1[0], position, b1[1], b1[6], b1[4], b2[0], b2[6], b2[4], str(same_length), synonymy
                outf.write("\t".join(line) + "\n")
