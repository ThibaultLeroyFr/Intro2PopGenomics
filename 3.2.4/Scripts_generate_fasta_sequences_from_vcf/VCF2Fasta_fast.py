#! /usr/bin/env python
from __future__ import print_function
import sys
import numpy as np
import optparse

# Function to find information in position in FORMAT field 
def findPosInFORMATStr(infoStr, str2ID):
	RETURNPOS = 0
	info = infoStr.split(":")
	for i in info:
		if i == str2ID:
			break
		RETURNPOS = RETURNPOS+1
	return(RETURNPOS)
	

parser = optparse.OptionParser()

parser.add_option('-q', action="store", dest="QUALThreshold", type="float")
parser.add_option('-m', '--min_cov',action="store", dest="minCoverageThreshold", type="int")
parser.add_option('-M', '--max_cov',action="store", dest="maxCoverageThreshold", type="int")
parser.add_option('-f', '--vcf_file',action="store", dest="namefp", type="string")

(options, args) = parser.parse_args()

minCoverageThreshold = int(options.minCoverageThreshold)
maxCoverageThreshold = int(options.maxCoverageThreshold)

##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
GQThreshold = float(options.QUALThreshold)

namefile = str(options.namefp)
fp = open(namefile)

############## MAIN ############
nameSeq = []
nameInd = []
sequenceDNA = ""


seqName = ""
countSNP = 0
countSite = 0
NumberOfSeq = ""
sites = []

pos = 0
oldpos = 0 

for line in fp:
	line = line.rstrip()
	
	# count number of individuals and store names in array
	if line[0:6] == "#CHROM" and len(nameInd) == 0:
		arrline = line.split()
		
		for i in arrline[9:]:
			indseq1 = i+".1"
			indseq2 = i+".2"
			nameSeq.append(indseq1)
			nameSeq.append(indseq2)
			nameInd.append(i)
		NumberOfSeq = len(nameSeq)
                print("There is "+str(NumberOfSeq/2)+" individuals")			
	if line[0] != "#":
		
		arrline = line.split()
		chromosome = arrline[0]
		pos = int(arrline[1])
		REF = str(arrline[3])
		ALT = str(arrline[4])
		site = []

		if (pos % 50000) == 0:
			print("Positions "+str(pos)+" of "+chromosome)

		
		if chromosome != seqName and seqName != "": # print seq when new chromosome found
			print("Print  " + seqName)
			fasta = open(seqName+".fst", "w")
			sites = np.ravel(sites)
			sequenceDNA = np.reshape(sites, newshape=(len(sites)/NumberOfSeq, NumberOfSeq))
			for i, seqName in enumerate(nameSeq):
					print(">"+seqName, file=fasta)
					print(''.join(map(str, sequenceDNA[1:,i])), file=fasta)
					
			seqName = chromosome
			sites = []
			
		if seqName == "":
			seqName = chromosome
		
		
		### MODIF PROPOSER PAR EMERIC START###############################
		## This was at the end before
		if oldpos+1 < pos: # if current site is not exactly one position after oldpos site put "N"
			numberSite = int(pos-oldpos)
			for j in range(1,numberSite):
				site = []
				for i in range(0, len(nameSeq)):
					site.append("N")
				#sequenceDNA = np.vstack((sequenceDNA, site))
                                sites.append(site)
                        site = []
		### MODIF PROPOSER PAR EMERIC STOP###############################
		
		
		if len(REF) > 1 or len(ALT) > 1 or ALT == "*" :  # exlcude indels
			for i in range(0, len(nameSeq)):
				site.append("N")
				
		else:

			covPos = findPosInFORMATStr(arrline[8], "DP") # find coverage position in INFO field

			if ALT == ".": # if no ALT allele ( = monorphic site)	
				if arrline[6] == "LowQual" : # exclude LowQual for HaplotypeCaller and GenotypeGVCFs SNP
					for i in range(0, len(nameSeq)):
						site.append("N")
				else:
					for i in range(0, len(nameInd)):
						
						ind = arrline[i+9]
						arrInd = ind.split(":")
						#print(str(len(arrInd)) + " " + str(covPos))
						if len(arrInd) < int(covPos+1): # if individual doesn't have DP field
							#print("here")
							site.append("N")
							site.append("N")
						else:
							cov = arrInd[covPos]
							if cov == ".":
								site.append("N")
								site.append("N")
								continue
							cov = int(cov)
							if cov > minCoverageThreshold and cov < maxCoverageThreshold:
								site.append(REF)
								site.append(REF)
							else:
								site.append("N")
								site.append("N")
						
			else:
				if arrline[6] != "PASS" : # exclude lowQ SNP
					for i in range(0, len(nameSeq)):
						site.append("N")
						
					#print("lowQ  " + line)
				else:
					#print("PASS  " + line)
					adPos = findPosInFORMATStr(arrline[8], "AD") #FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
					gqPos = findPosInFORMATStr(arrline[8], "GQ") #FORMAT=<ID=FT,Number=.,Type=String,Description="Genotype-level filter">
					for i in range(0, len(nameInd)):
						ind = arrline[i+9]
						arrInd = ind.split(":")
						if len(arrInd) < int(gqPos+1): # if individual doesn't have GQ field
							site.append("N")
							site.append("N")
							continue
							
						GQ = arrInd[gqPos]
						if GQ == "." :
							site.append("N")
							site.append("N")
							continue
						arrInd = ind.split(":")
						ad = arrInd[adPos].split(",")
						if int(ad[0]) > 0:
							site.append(REF)
						else:
							site.append(ALT)
						if int(ad[1]) > 0:
							site.append(ALT)
						else:
							site.append(REF)
						
		


                sites.append(site)
                
                ## Test if the site added are of correct length
                #sit = np.ravel(sites)
                #if not (float(len(sit))/float(NumberOfSeq)).is_integer()  :
                #        print(chromosome)
                #        print("Site "+str(len(site)))
                #        print("Sites "+str(len(sit)))
                #        print(str(pos))
                #        sys.exit("NOOOO :-)")

                
		oldpos = pos
		
		

print("Print  " + seqName)
fasta = open(seqName+".fst", "w")
sites = np.ravel(sites)
sequenceDNA = np.reshape(sites, newshape=(len(sites)/NumberOfSeq, NumberOfSeq))
for i, seqName in enumerate(nameSeq):
		print(">"+seqName, file=fasta)
		print(''.join(map(str, sequenceDNA[1:,i])), file=fasta)


#print "Number of SNP "+str(count)

fp.close()

