#!/usr/bin/env python
# Modified by TL 012217 

import sys
import numpy as np

## Compute the depth of coverage per individuals.

fp = open(sys.argv[1])

nameInd = []
depthOfCoverage = ""

for line in fp:
	line = line.rstrip()
	
	# count number of individuals and store names in array
	if line[0:6] == "#CHROM":
		arrline = line.split()
		for i in arrline[9:]:
			nameInd.append(i)
		depthOfCoverage = nameInd
	if line[0] != "#":
		arrline = line.split()
		chromosome = arrline[0]
		pos = int(arrline[1])
		# loop other the individuals
                if arrline[6] != "PASS":
                    #print line
		    cov = []
                    qual = []
		    covPos = 0
                    qualPos = 0
		    info = arrline[8].split(":")
		    for i in info:
			    if i == "DP":
                                break
                            covPos= covPos+1
                    for i in info:
                            if i == "GQ":
                                break
                            qualPos = qualPos+1

                    #print (covPos,qualPos)
                    #print chromosome,pos,info
		    for i in range(0, len(nameInd)):
                            ID = nameInd[i]
			    ind = arrline[i+9]
			    arrInd = ind.split(":")
                            if len(arrInd) > qualPos: # virer les cas ou les ind n ont pas assez de champs
                                mycov = arrInd[covPos]
                                myqual = arrInd[qualPos]
                                myline = str(ID) + "\t" + str(chromosome) + "\t" + str(pos) + "\t" + str(mycov) + "\t" + str(myqual)
		                #cov.append(myline) # add quality
                                print(myline)
			
		# store coverage in a matrix
		#depthOfCoverage = np.vstack((depthOfCoverage, cov))


#for i, element in enumerate(depthOfCoverage):
#	print '\n'.join(element) # \t initialement

fp.close()

