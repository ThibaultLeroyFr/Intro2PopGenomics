# TL - 260319
# this script generates maf allele freq from all SNPs in all individuals

import sys
import re

file=sys.argv[1]
fileop=open(file)

filelistop=open(file+".list")

out=open(file+".Og.maf","w")

dico = {}

for line in fileop.readlines(): # geno file [with col1=indID col2=scaffID col3=pos ...]
    line=line.strip()
    splittedline=line.split("\t")
    individual=splittedline[0]
    scaffID=splittedline[1]
    pos=splittedline[2]
    #print(scaffID,pos)
    key = (individual + scaffID + pos)
    dico[key] = line

for line2 in filelistop.readlines(): # list of files
    total=""
    line2=line2.strip()
    splittedline2=line2.split("\t")
    scaffID2=splittedline2[0]
    pos2=splittedline2[1]
    pattern24=("Og_ERR2008921" + scaffID2 + pos2)
    pattern25=("Og_ERR2008922" + scaffID2 + pos2)
    pattern26=("Og_ERR2008923" + scaffID2 + pos2)
    pattern27=("Og_ERR2008927" + scaffID2 + pos2)
    pattern28=("Og_ERR2008934" + scaffID2 + pos2)
    pattern29=("Og_ERR2008935" + scaffID2 + pos2)
    pattern30=("Og_ERR2008946" + scaffID2 + pos2)
    pattern31=("Og_ERR2008953" + scaffID2 + pos2)
    pattern32=("Og_ERR2008955" + scaffID2 + pos2)
    pattern33=("Og_ERR2008959" + scaffID2 + pos2)
    pattern34=("Og_ERR2008963" + scaffID2 + pos2)
    pattern35=("Og_ERR2008974" + scaffID2 + pos2)
    pattern36=("Og_ERR2008980" + scaffID2 + pos2)
    pattern37=("Og_ERR2008984" + scaffID2 + pos2)
    pattern38=("Og_ERR2008998" + scaffID2 + pos2)
    pattern39=("Og_ERR2008999" + scaffID2 + pos2)
    pattern40=("Og_ERR2009002" + scaffID2 + pos2)
    pattern41=("Og_ERR2009006" + scaffID2 + pos2)
    pattern42=("Og_ERR2009023" + scaffID2 + pos2)
    pattern43=("Og_ERR2009026" + scaffID2 + pos2)
    pattern44=("Og_ERR2009027" + scaffID2 + pos2)
    pattern45=("Og_ERR2009041" + scaffID2 + pos2)
    pattern46=("Og_ERR2009051" + scaffID2 + pos2)
    pattern47=("Og_ERR2009065" + scaffID2 + pos2)
    pattern48=("Og_ERR2009075" + scaffID2 + pos2)
    if (dico.has_key(pattern24)):
        #print dico[pattern4]
        total=total+dico[pattern24]
    if (dico.has_key(pattern25)):
        #print dico[pattern5]
        total=total+dico[pattern25]
    if (dico.has_key(pattern26)):
        #print dico[pattern5]
        total=total+dico[pattern26]
    if (dico.has_key(pattern27)):
        #print dico[pattern5]
        total=total+dico[pattern27]
    if (dico.has_key(pattern28)):
        #print dico[pattern1]
        total=total+dico[pattern28]
    if (dico.has_key(pattern29)):
        #print dico[pattern2]
        total=total+dico[pattern29]
    if (dico.has_key(pattern30)):
        #print dico[pattern3]
        total=total+dico[pattern30]
    if (dico.has_key(pattern31)):
        #print dico[pattern1]
        total=total+dico[pattern31]
    if (dico.has_key(pattern32)):
        #print dico[pattern2]
        total=total+dico[pattern32]
    if (dico.has_key(pattern33)):
        #print dico[pattern3]
        total=total+dico[pattern33]
    if (dico.has_key(pattern34)):
        #print dico[pattern4]
        total=total+dico[pattern34]
    if (dico.has_key(pattern35)):
        #print dico[pattern5]
        total=total+dico[pattern35]
    if (dico.has_key(pattern36)):
        #print dico[pattern5]
        total=total+dico[pattern36]
    if (dico.has_key(pattern37)):
        #print dico[pattern5]
        total=total+dico[pattern37]
    if (dico.has_key(pattern38)):
        #print dico[pattern1]
        total=total+dico[pattern38]
    if (dico.has_key(pattern39)):
        #print dico[pattern2]
        total=total+dico[pattern39]
    if (dico.has_key(pattern40)):
        #print dico[pattern3]
        total=total+dico[pattern40]
    if (dico.has_key(pattern41)):
        #print dico[pattern1]
        total=total+dico[pattern41]
    if (dico.has_key(pattern42)):
        #print dico[pattern2]
        total=total+dico[pattern42]
    if (dico.has_key(pattern43)):
        #print dico[pattern3]
        total=total+dico[pattern43]
    if (dico.has_key(pattern44)):
        #print dico[pattern4]
        total=total+dico[pattern44]
    if (dico.has_key(pattern45)):
        #print dico[pattern5]
        total=total+dico[pattern45]
    if (dico.has_key(pattern46)):
        #print dico[pattern5]
        total=total+dico[pattern46]
    if (dico.has_key(pattern47)):
        #print dico[pattern5]
        total=total+dico[pattern47]
    if (dico.has_key(pattern48)):
        #print dico[pattern1]
        total=total+dico[pattern48]
    counthomoref=0
    counthomoalt=0
    counthetero=0
    countmissing=0
    allfreqalt=float(0)
    maf=float(0)
    #pattern=str(0/0)
    counthomoref=sum(1 for match in re.finditer(r"0/0", total))
    counthomoalt=sum(1 for match in re.finditer(r"1/1", total))
    counthetero=sum(1 for match in re.finditer(r"0/1", total))
    countmissing=sum(1 for match in re.finditer(r"\./\.", total))
    if countmissing == 0:
        allfreqalt=float((2*float(counthomoalt)+int(counthetero))/(2*(int(counthomoref) + int(counthomoalt) + int(counthetero))))
        if allfreqalt <= 0.5:
            maf=allfreqalt
        else:
            maf=1-allfreqalt
    else:
        allfreqalt="NA"
        maf="NA"
    #for pattern in total:
    #    counthomoref+=1
    myline=str(scaffID2)+"\t"+str(pos2)+"\t"+str(counthomoref)+"\t"+str(counthomoalt)+"\t"+str(counthetero)+"\t"+str(countmissing)+"\t"+str(allfreqalt)+"\t"+str(maf)+"\n"
    #print(myline)
    out.write(myline)

out.write.close()
    


