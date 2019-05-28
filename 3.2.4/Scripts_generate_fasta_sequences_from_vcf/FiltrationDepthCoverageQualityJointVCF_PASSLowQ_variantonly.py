#!/usr/bin/env python
# Modified by TL 012217 

import sys
import numpy as np
import re

## Extract champs de qualite joinVCF
# pour lheure je me base sur des fichiers generes pour les variants uniquements avec ce genre de ligne: grep -v "#" joint_Fsemitorquata_bwa_mem_mdup_filtered.vcf | grep -v "LowQual" | awk '$7 != "." {print $1"        "$2"    "$3"    "$4"    "$5"    "$6"    "$8"    "$7}' > joint_Fsemitorquata_bwa_mem_mdup_filtered.vcf.qual.withPASS

fp = open(sys.argv[1])

nameInd = []
depthOfCoverage = ""

print("scaffold	pos	type	ref	alt	label   QUAL	AC	AF	AN	DP	ExcessHet	MQ	QD")

for line in fp:
	line = line.rstrip()
	# ignore line with # if any
	if line[0] != "#":
		arrline = line.split()
		chromosome = arrline[0]
		pos = int(arrline[1])
                ref = arrline[3]
                alt = arrline[4]
                qual=arrline[5]
                label=arrline[7] # PASS or lowQ
		# loop other the individuals
                testMQ = 0
                testQD = 0
		myAC = []
                myAP = []
                myAN = []
                myDP = []
                myExcess = []
                myMQ = []
                myQD = []
		ACPos = 0
                AFPos = 0
                ANPos = 0
                DPPos = 0
                ExcessPos = 0
                MQPos = 0
                QDPos = 0
                #info = arrline[6].split(";").split("=")
                info = re.split(r'[;=]', arrline[6])
                #print info
		for i in info:
			if i == "AC":
                          break
                        ACPos = ACPos+1
                for i in info:
                        if i == "AF":
                            break
                        AFPos = AFPos+1
                for i in info:
                        if i == "AN":
                            break
                        ANPos = ANPos+1
                for i in info:
                        if i == "DP":
                            break
                        DPPos = DPPos+1
                for i in info:
                        if i == "ExcessHet":
                            break
                        ExcessPos = ExcessPos+1
                for i in info:
                        if i == "MQ":
                            testMQ = 1
                            break
                        MQPos = MQPos+1
                for i in info:
                        if i == "QD":
                            testQD = 1
                            break
                        QDPos = QDPos+1
                # ADD 1
                ACPos = ACPos+1
                AFPos = AFPos+1
                ANPos = ANPos+1
                DPPos = DPPos+1
                ExcessPos = ExcessPos+1
                MQPos = MQPos+1
                QDPos = QDPos+1
                #print (covPos,qualPos)
                #print chromosome,pos,info
		myAC = info[ACPos]
                myAF = info[AFPos]
                myAN = info[ANPos]
                myDP = info[DPPos]
                myExcess = info[ExcessPos]
                if testMQ == 1:
                    myMQ = info[MQPos]
                else:
                    myMQ = "NA"
                if testQD == 1: 
                    myQD = info[QDPos]
                else:
                    myQD = "NA"
                if len(ref) > 1 or len(alt) > 1: # if more than two characters in columns ref and alt
                    if "," in ref or "," in alt: # test if more than 2 alleles or several nucleotides (INDEL)
                        myline = str(chromosome) + "\t" + str(pos) + "\tMULTIPLE\t" + str(ref) + "\t" + str(alt) + "\t" + str(label) + "\t" + str(qual) + "\t" + str(myAC) + "\t" + str(myAF) + "\t" + str(myAN) +"\t"+ str(myDP) + "\t" + str(myExcess)+"\t"+str(myMQ)+"\t"+str(myQD)
                    else:
                        myline = str(chromosome) + "\t" + str(pos) + "\tINDEL\t" + str(ref) + "\t" + str(alt) + "\t" + str(label) + "\t" + str(qual) + "\t"+ str(myAC) + "\t" + str(myAF) + "\t" + str(myAN) +"\t"+ str(myDP) + "\t" + str(myExcess)+"\t"+str(myMQ)+"\t"+str(myQD)
                else:
                    myline = str(chromosome) + "\t" + str(pos) + "\tSNP\t" + str(ref) + "\t" + str(alt) + "\t" + str(label)+ "\t" +  str(qual) + "\t" + str(myAC) + "\t" + str(myAF) + "\t" + str(myAN) +"\t"+ str(myDP) + "\t" + str(myExcess)+"\t"+str(myMQ)+"\t"+str(myQD)
		#cov.append(myline) # add quality
                print(myline)
			
		# store coverage in a matrix
		#depthOfCoverage = np.vstack((depthOfCoverage, cov))


#for i, element in enumerate(depthOfCoverage):
#	print '\n'.join(element) # \t initialement

fp.close()

