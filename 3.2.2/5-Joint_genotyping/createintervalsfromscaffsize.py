#!/usr/bin/env python
#$Id$


### usage: python createintervalsfromscaffsize.py [FASTA_FILE] [NB_INTERVALS]
### eg: python createintervalsfromscaffsize.py ZoBo_SOAP_SSPACElongread_PE_MP_PacBio.fasta.scaffsize 100


import os
import re
import string
import sys
import glob

# rm file if exists
try:
	filelist=glob.glob("interval*.txt")
	for file in filelist:
		os.remove(file)
except OSError:
	pass

# first round, total scaff size
file1 = sys.argv[1] # open the s
intervals = sys.argv[2] # open the s
file1_stream = open(file1)
totalsize=0
for line1 in file1_stream.readlines(): 
	line1 = line1.replace('\n','')
	splitted_line1 = line1.split('\t')
	scaffID = splitted_line1[0] 
	size = int(splitted_line1[1]) 
	totalsize += size

print "assembly=",totalsize
print "creating ",intervals," intervals"


limit=totalsize / int(intervals)
print "limit for intervals",limit

file1_stream = open(file1)

# second round, create intervals
countfile=0
countinterval=1
currentinterval="interval_"+str(countinterval)
previoustotalsize=0
positionremaining=0
previousID=""
previouslimit=1
currentlimit=limit
remaining=0
totalremaining=0
loopwhile=0
for line1 in file1_stream.readlines(): 
	line1 = line1.replace('\n','')
	splitted_line1 = line1.split('\t')
	scaffID = splitted_line1[0] 
	size = int(splitted_line1[1]) 
	lastfirstposition=previoustotalsize - positionremaining + 1
	currentlimit = limit - positionremaining # next step will be the total limit for the interval minus the number of bases already use to end the previous scaff
	if loopwhile == 1 and totalremaining > 0:
		previouslimit=1
		#file2.close()
		countfile+=1
		myfile="scatter"+str(countfile)+".intervals"
		file2=open(myfile,"w")
                toprint=previousID+":"+str(lastfirstposition)+"-"+str(previoustotalsize)+"\n"# print the rest of the line
		#print previousID,"\t",str(lastfirstposition),"\t",str(previoustotalsize),"interval_",countinterval,"\tfile",countfile,"\ttotalremaining=",totalremaining
		#totalremaining=0
		countinterval+=1
		currentinterval="interval_"+str(countinterval)
		file2.write(toprint)
	if size > currentlimit : # if scaffold greater the limit # previously and size > 100000
		while size > currentlimit: # loop to cut the scaffolds into intervals 
			if previouslimit == 1 and totalremaining > 0:
				myfile="scatter"+str(countfile)+".intervals"
				file2=open(myfile,"a")
				toprint=scaffID+":"+str(previouslimit)+"-"+str(currentlimit)+"\n"
				#print scaffID,"\t",str(previouslimit),"\t",str(currentlimit),"interval_",countinterval,"\tfile",countfile,"\ttotalremaining=",totalremaining
				countinterval+=1
				currentinterval="interval_"+str(countinterval)
				file2.write(toprint)
				file2.close()
				previouslimit=currentlimit+1
				currentlimit+=limit
				#totalremaining=0
			else:
				countfile+=1
				myfile="scatter"+str(countfile)+".intervals"
				file2=open(myfile,"a")
                                toprint=scaffID+":"+str(previouslimit)+"-"+str(currentlimit)+"\n"
				#print scaffID,"\t",str(previouslimit),"\t",str(currentlimit),"interval_",countinterval,"\tfile",countfile,"\ttotalremaining=",totalremaining
				countinterval+=1
				currentinterval="interval_"+str(countinterval)
				file2.write(toprint)
				previouslimit+=limit
				currentlimit+=limit
				totalremaining=0
				file2.close()
			loopwhile=1
		previoustotalsize=size
		positionremaining=size-(currentlimit-limit) # remaining
		totalremaining+=positionremaining
	else: # for all other scaffolds, add the length of the current scaffold and release only if the total length is greater than the limit for each interval
		if countfile == 0: # if the first scaffold is lower than the limit, open a file
			countfile+=1
			myfile="scatter"+str(countfile)+".intervals"
			file2=open(myfile,"w")
		loopwhile=0
		totalremaining+=size
		if totalremaining > limit:
                        toprint=scaffID+":1-"+str(size)+"\n"
			file2.write(toprint)
			#print scaffID,"\t1\t",str(size),"interval_",countinterval,"\tfile",countfile,"\tremaining=",remaining,"\ttotalremaining=",totalremaining
			file2.close()
			totalremaining=0
			countfile+=1
			countinterval+=1
			currentinterval="interval_"+str(countinterval)
			myfile="scatter"+str(countfile)+".intervals"
			file2=open(myfile,"w")
		else:
                        toprint=scaffID+":1-"+str(size)+"\n"
			file2.write(toprint)
			#print scaffID,"\t1\t",str(size),"interval_",countinterval,"\tfile",countfile,"\tremaining=",remaining,"\ttotalremaining=",totalremaining
			countinterval+=1
			currentinterval="interval_"+str(countinterval)
	
	#print positionremaining
	previousID=scaffID

file2.close()

	
	
