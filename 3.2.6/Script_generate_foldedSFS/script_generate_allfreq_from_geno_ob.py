# TL - 260319
# this script generates maf allele freq from all SNPs in all individuals

import sys
import re

file=sys.argv[1]
fileop=open(file)

filelistop=open(file+".list")

out=open(file+".Ob.maf","w")

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
    pattern1=("Ob_ERR2008849" + scaffID2 + pos2)
    pattern2=("Ob_ERR2008850" + scaffID2 + pos2)
    pattern3=("Ob_ERR2008851" + scaffID2 + pos2)
    pattern4=("Ob_ERR2008852" + scaffID2 + pos2)
    pattern5=("Ob_ERR2008853" + scaffID2 + pos2)
    pattern6=("Ob_ERR2008855" + scaffID2 + pos2)
    pattern7=("Ob_ERR2008856" + scaffID2 + pos2)
    pattern8=("Ob_ERR2008857" + scaffID2 + pos2)
    pattern9=("Ob_ERR2008858" + scaffID2 + pos2)
    pattern10=("Ob_ERR2008859" + scaffID2 + pos2)
    pattern11=("Ob_ERR2008860" + scaffID2 + pos2)
    pattern12=("Ob_ERR2008861" + scaffID2 + pos2)
    pattern13=("Ob_ERR2008862" + scaffID2 + pos2)
    pattern14=("Ob_ERR2008863" + scaffID2 + pos2)
    pattern15=("Ob_ERR2008884" + scaffID2 + pos2)
    pattern16=("Ob_ERR2008894" + scaffID2 + pos2)
    pattern17=("Ob_ERR2008910" + scaffID2 + pos2)
    pattern18=("Ob_ERR2008914" + scaffID2 + pos2)
    pattern19=("Ob_ERR2009084" + scaffID2 + pos2)
    pattern20=("Ob_ERR2009085" + scaffID2 + pos2)
    pattern21=("Ob_ERR2009086" + scaffID2 + pos2)
    pattern22=("Ob_ERR2009090" + scaffID2 + pos2)
    pattern23=("Ob_ERR2009096" + scaffID2 + pos2)
    if (dico.has_key(pattern1)):
        #print dico[pattern1]
        total=total+dico[pattern1]
    if (dico.has_key(pattern2)):
        #print dico[pattern2]
        total=total+dico[pattern2]
    if (dico.has_key(pattern3)):
        #print dico[pattern3]
        total=total+dico[pattern3]
    if (dico.has_key(pattern4)):
        #print dico[pattern4]
        total=total+dico[pattern4]
    if (dico.has_key(pattern5)):
        #print dico[pattern5]
        total=total+dico[pattern5]
    if (dico.has_key(pattern6)):
        #print dico[pattern5]
        total=total+dico[pattern6]
    if (dico.has_key(pattern7)):
        #print dico[pattern5]
        total=total+dico[pattern7]
    if (dico.has_key(pattern8)):
        #print dico[pattern1]
        total=total+dico[pattern8]
    if (dico.has_key(pattern9)):
        #print dico[pattern2]
        total=total+dico[pattern9]
    if (dico.has_key(pattern10)):
        #print dico[pattern3]
        total=total+dico[pattern10]
    if (dico.has_key(pattern11)):
        #print dico[pattern1]
        total=total+dico[pattern11]
    if (dico.has_key(pattern12)):
        #print dico[pattern2]
        total=total+dico[pattern12]
    if (dico.has_key(pattern13)):
        #print dico[pattern3]
        total=total+dico[pattern13]
    if (dico.has_key(pattern14)):
        #print dico[pattern4]
        total=total+dico[pattern14]
    if (dico.has_key(pattern15)):
        #print dico[pattern5]
        total=total+dico[pattern15]
    if (dico.has_key(pattern16)):
        #print dico[pattern5]
        total=total+dico[pattern16]
    if (dico.has_key(pattern17)):
        #print dico[pattern5]
        total=total+dico[pattern17]
    if (dico.has_key(pattern18)):
        #print dico[pattern1]
        total=total+dico[pattern18]
    if (dico.has_key(pattern19)):
        #print dico[pattern2]
        total=total+dico[pattern19]
    if (dico.has_key(pattern20)):
        #print dico[pattern3]
        total=total+dico[pattern20]
    if (dico.has_key(pattern21)):
        #print dico[pattern1]
        total=total+dico[pattern21]
    if (dico.has_key(pattern22)):
        #print dico[pattern2]
        total=total+dico[pattern22]
    if (dico.has_key(pattern23)):
        #print dico[pattern3]
        total=total+dico[pattern23]
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
    


