#! /usr/bin/env python
from __future__ import print_function
import sys
import numpy as np
import optparse

## Filter variant from a GVCF file obtained by genotyping one individual (not from joint genotyping using "GenotypeGVCFs" programme)


# Function to find information in position in FORMAT field 
def findPosInFORMATStr(infoStr, str2ID):
	RETURNPOS = 0
	info = infoStr.split(":")
	for i in info:
		if i == str2ID:
			break
		RETURNPOS = RETURNPOS+1
	return(RETURNPOS)
	
	
# Function to find information in INFO field  and return the values
def returnValueInINFOStr(infoStr, str2ID):
	RETURNVALUE = -999
	info = infoStr.split(";")
	for i in info:
		arr = i.split("=")
		if str2ID == str(arr[0]):
			RETURNVALUE = float(arr[1])
	return(RETURNVALUE)
	
	
############# OPTION ###################################################
	
parser = optparse.OptionParser()
parser.add_option('-q', action="store", dest="QUALThreshold", help="QUALThreshold; QD; if QD < X -> lowQ", type="float")
parser.add_option('-s', '--FS', action="store", dest="FS", help="FS; FS > X  -> lowQ; Phred-scaled p-value using Fisher's exact test to detect strand bias", type="float")
parser.add_option('-m', '--MQ', action="store", dest="MQ", help="MQ; MQ < X -> lowQ; RMS Mapping Quality", type="float")
parser.add_option('-n', '--MQRankSum', action="store", dest="MQRankSum", help="MQRankSum; MQRankSum < X -> lowQ; Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities", type="float")
parser.add_option('-r', '--ReadPosRankSum', action="store", dest="ReadPosRankSum", help="ReadPosRankSum < X -> lowQ; ReadPosRankSum, Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias", type="float")
parser.add_option('-w', '--RAW_MQ', action="store", dest="RAW_MQ", help="RAW_MQ; RAW_MQ < X -> lowQ", type="float")
parser.add_option('-f', '--vcf_file',action="store", dest="namefp", type="string")

(options, args) = parser.parse_args()

QD_Thres = float(options.QUALThreshold)
FS_Thres = float(options.FS)
MQ_Thres = float(options.MQ)
MQRankSum_Thres = float(options.MQRankSum)
RAW_MQ_Thres = float(options.RAW_MQ)
ReadPosRankSum_Thres = float(options.ReadPosRankSum)

namefile = str(options.namefp)
fp = open(namefile)


for line in fp:
	line = line.rstrip()
	
	# count number of individuals and store names in array
	if line[0] == "#":
		print(line)
		continue
			
	if line[0] != "#":
		arrline = line.split()
		chromosome = arrline[0]
		pos = int(arrline[1])
		REF = str(arrline[3])
		ALT = str(arrline[4])	
		
		
		if ALT == "<NON_REF>" or ALT == ".": # if molymorphic site
			print(line)
			continue
			
		ALT = ALT.replace(",<NON_REF>", "NN")
                #ALT = ALT.replace(",.", "")
		if len(REF) > 1 or len(ALT) > 1 or ALT == "*" :  # exlcude indels 
			print(line)
			continue
			
		#print("HERE")
		PASS = "PASS"
		QD = returnValueInINFOStr(arrline[7], "QD")
		if QD < QD_Thres and QD != -999:
			PASS = "lowQ"
		#print("QD "+str(QD))
		
		FS = returnValueInINFOStr(arrline[7], "FS")
		if FS > FS_Thres and QD != -999:
			PASS = "lowQ"
		#print("FS "+str(FS))
		
		MQ = returnValueInINFOStr(arrline[7], "MQ")
                if MQ < MQ_Thres and MQ != -999:
                        PASS = "lowQ"
                elif MQ == "NaN":
			PASS = "lowQ"
		#print("MQ "+str(MQ))
		if MQ == -999: # if MQ is absent use RAW_MQ
			RAW_MQ = returnValueInINFOStr(arrline[7], "RAW_MQ")
			if RAW_MQ < RAW_MQ_Thres:
				PASS = "lowQ"
				
		MQRankSum = returnValueInINFOStr(arrline[7], "MQRankSum")
		if MQRankSum < MQRankSum_Thres and MQRankSum != -999:
			PASS = "lowQ"
		#print("MQRankSum "+str(MQRankSum))
		
		ReadPosRankSum = returnValueInINFOStr(arrline[7], "ReadPosRankSum")
		if ReadPosRankSum < ReadPosRankSum_Thres and ReadPosRankSum != -999:
			PASS = "lowQ"
		#print("ReadPosRankSum "+str(ReadPosRankSum))
		#print(arrline[7])
		
		c = 0
		newline =""
		for it in arrline:
			c = c + 1
			if c == 1: 
				newline = str(it)
				continue
			if c == 7:
				newline = newline + "\t" + PASS
				continue
				
			newline = newline + "\t" + str(it)
		print(newline )
fp.close()

