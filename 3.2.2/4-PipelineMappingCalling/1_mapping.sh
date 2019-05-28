#! /usr/bin/env bash

######### 

## define global variables
acc=$1
readsDIR=$2
nfq1=$3
nfq2=$4
nfqu1=$5
nfqu2=$6
refgen=$7
nthr=$8

fq1=${readsDIR}/${nfq1}
fq2=${readsDIR}/${nfq2}
fqu1=${readsDIR}/${nfqu1}
fqu2=${readsDIR}/${nfqu2}

tmpdir=$(echo "/tmp/")

### set the file path to bwa
bwa=$(echo "/sandbox/users/tleroy/Software/bwa-0.7.17/bwa") # bwa version we used: 0.7.17
picard=$(echo "/sandbox/users/tleroy/Software/picard-tools-1.140/picard.jar") # picard version we used: 1.140 (requires java )
samtools=$(echo "/sandbox/users/tleroy/Software/samtools-1.3.1/samtools") # samtools version we used: 1.3.1


################################# No changes expected below #################################
####### MAPPING #########
## OPTION BWA
# option -R = reads group

rgstring="@RG\tID:"$acc"\tDT:2017-01-25\tLB:lib-$acc\tPL:ILLUMINA\tSM:$acc"


bwape=${acc}_bwa_mem
bwase1=${acc}_bwa_mem_sing1
bwase2=${acc}_bwa_mem_sing2


cmd="$bwa mem -R \"${rgstring}\" \
	-M \
	-t ${nthr} \
	${refgen} \
	${fq1} ${fq2} > ${bwape}.sam"
 
# display the full command with expanded variables
echo "# the command will be:"
echo ${cmd}
 
# execute the command
eval ${cmd} 2>${acc}_bwa_mem.err

# bwa se unpaired 1
cmd="$bwa mem -R \"${rgstring}\" \
	-M \
	-t ${nthr} \
	${refgen} \
	${fqu1}  > ${bwase1}.sam"
# display the full command with expanded variables
echo "# the command will be:"
echo ${cmd}
# execute the command
eval ${cmd} 2>${acc}_bwa_mem_sing1.err

# bwa se unpaired 2
cmd="$bwa mem -R \"${rgstring}\" \
	-M \
	-t ${nthr} \
	${refgen} \
	${fqu2}  > ${bwase2}.sam"
# display the full command with expanded variables
echo "# the command will be:"
echo ${cmd}
# execute the command
eval ${cmd} 2>${acc}_bwa_mem_sing2.err


################# FILTERING  ##############

#Remove unmapped and filter on quality
echo "Remove unmapped and filter on quality"
$samtools view -S -q 20  ${bwape}.sam > ${bwape}.filtered.sam
grep '^@' ${bwape}.sam > ${bwape}_headers_sam
cat ${bwape}_headers_sam ${bwape}.filtered.sam > ${bwape}.sam
rm ${bwape}_headers_sam ${bwape}.filtered.sam

$samtools view -S -q 20  ${bwase1}.sam > ${bwase1}.filtered.sam
grep '^@' ${bwase1}.sam > ${bwase1}_headers_sam
cat ${bwase1}_headers_sam ${bwase1}.filtered.sam > ${bwase1}.sam
rm ${bwase1}_headers_sam ${bwase1}.filtered.sam

$samtools view -S -q 20 ${bwase2}.sam > ${bwase2}.filtered.sam
grep '^@' ${bwase2}.sam > ${bwase2}_headers_sam
cat ${bwase2}_headers_sam ${bwase2}.filtered.sam > ${bwase2}.sam
rm ${bwase2}_headers_sam ${bwase2}.filtered.sam


 
##### convert to sorted bam, index & cleanup
echo "convert to sorted bam, index & cleanup"
java -Xmx4g -Djava.io.tmpdir=$tmpdir -jar $picard SortSam I=${bwape}.sam O=${bwape}.bam SO=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmpdir
rm ${bwape}.sam
java -Xmx4g -Djava.io.tmpdir=$tmpdir -jar $picard SortSam I=${bwase1}.sam O=${bwase1}.bam SO=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmpdir
rm ${bwase1}.sam
java -Xmx4g -Djava.io.tmpdir=$tmpdir -jar $picard SortSam I=${bwase2}.sam O=${bwase2}.bam SO=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmpdir
rm ${bwase2}.sam


$samtools merge ${bwape}_merged.bam ${bwape}.bam ${bwase1}.bam ${bwase2}.bam
ls -lh ${bwape}.bam >>${bwape}.log
mv ${bwape}_merged.bam ${bwape}.bam 
rm ${bwase1}.bam ${bwase2}.bam
ls -lh ${bwape}.bam >>${bwape}.log

##### mark duplicates with Picard & index
# READ_NAME_REGEX=null # read name does not contain information physical location of the cluster in the lane
echo "mark duplicates with Picard & index"
java -Xmx4g -Djava.io.tmpdir=$tmpdir -jar $picard MarkDuplicates \
   I=${bwape}.bam \
   O=${bwape}_mdup.bam \
   MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
   CREATE_INDEX=true \
   VALIDATION_STRINGENCY=LENIENT \
   READ_NAME_REGEX=null  \
   METRICS_FILE=${bwape}_duplicate_metrics.txt \
   TMP_DIR=$tmpdir
   
rm ${bwape}.bam
rm ${bwape}.bai
   
#java -jar /bigvol/benoit/software/picard-tools-1.140/picard.jar BuildBamIndex \ 
#   INPUT=${bwape}_mdup.bam
	
##### get basic stats from the resulting bam file with Samtools
echo "Get basic stats from the resulting bam file with Samtools"
$samtools flagstat ${bwape}_mdup.bam >${bwape}_mdup_flagstat.txt
	

echo "# finished mapping ${pfx} reads with BWA mem"
