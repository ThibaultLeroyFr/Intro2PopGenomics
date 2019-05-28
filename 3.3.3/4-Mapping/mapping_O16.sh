# TL - 05/06/15 - thibault.leroy@pierroton.inra.fr
# pipeline mapping using bowtie2, picards & samtools
### qsub parameters
#$ -pe parallel_smp 8
#$ -m abe

#ERR3284937
#ERR3284938
#ERR3284939
#ERR3284940
#
##################### TO BE EDITED ###############################
ref=$(echo /work2/project/genoak/tleroy/ref/Qrob_PM1N) # set the path to your reference # bowtie2-build ref needed!
root=$(echo /work2/project/genoak/tleroy/SessileOak/O16)
root1=$(echo /work2/project/genoak/tleroy/SessileOak/ERR3284937) # set the path1
root2=$(echo /work2/project/genoak/tleroy/SessileOak/ERR3284938) # set the path2
root3=$(echo /work2/project/genoak/tleroy/SessileOak/ERR3284939) # set the path3
root4=$(echo /work2/project/genoak/tleroy/SessileOak/ERR3284940) # set the path4
############# NO CHANGES EXPECTED BELOW ##########################
# pair_end file name edition
change1=$(echo "$root1.1")
change2=$(echo "$root1.2")
change3=$(echo "$root2.1")
change4=$(echo "$root2.2")
change5=$(echo "$root3.1")
change6=$(echo "$root3.2")
change7=$(echo "$root4.1")
change8=$(echo "$root4.2")
pair_end1=$(echo "$change1""_clean.fastq.gz") # assuming paired trimmed reads: ERR3284937.1.clean.fastq.gz
pair_end2=$(echo "$change2""_clean.fastq.gz") # assuming paired trimmed reads: ERR3284937.2.clean.fastq.gz
pair_end3=$(echo "$change3""_clean.fastq.gz")
pair_end4=$(echo "$change4""_clean.fastq.gz")
pair_end5=$(echo "$change5""_clean.fastq.gz")
pair_end6=$(echo "$change6""_clean.fastq.gz")
pair_end7=$(echo "$change7""_clean.fastq.gz")
pair_end8=$(echo "$change8""_clean.fastq.gz")
# bowtie2 analysis, outputs bam files
/usr/local/bioinfo/src/bowtie/current2/bowtie2 -p8 -k1 -q --sensitive -x $ref -1 $pair_end1 -2 $pair_end2 -t --un-gz $root1.v2.3 | samtools view -Shb /dev/stdin > $root1.v2.3.bam
/usr/local/bioinfo/src/bowtie/current2/bowtie2 -p8 -k1 -q --sensitive -x $ref -1 $pair_end3 -2 $pair_end4 -t --un-gz $root2.v2.3 | samtools view -Shb /dev/stdin > $root2.v2.3.bam
/usr/local/bioinfo/src/bowtie/current2/bowtie2 -p8 -k1 -q --sensitive -x $ref -1 $pair_end5 -2 $pair_end6 -t --un-gz $root3.v2.3 | samtools view -Shb /dev/stdin > $root3.v2.3.bam
/usr/local/bioinfo/src/bowtie/current2/bowtie2 -p8 -k1 -q --sensitive -x $ref -1 $pair_end7 -2 $pair_end8 -t --un-gz $root4.v2.3 | samtools view -Shb /dev/stdin > $root4.v2.3.bam
# print unmapped reads
samtools view -f 4 -b $root1.v2.3.bam > $root1.v2.3.bam.unmapped
samtools view -f 4 -b $root2.v2.3.bam > $root2.v2.3.bam.unmapped
samtools view -f 4 -b $root3.v2.3.bam > $root3.v2.3.bam.unmapped
samtools view -f 4 -b $root4.v2.3.bam > $root4.v2.3.bam.unmapped
# picard sort
java -Xmx4g -jar /usr/local/bioinfo/src/picard-tools/current/SortSam.jar INPUT=$root1.v2.3.bam OUTPUT=$root1.v2.3.bam.pisorted SORT_ORDER=coordinate
java -Xmx4g -jar /usr/local/bioinfo/src/picard-tools/current/SortSam.jar INPUT=$root2.v2.3.bam OUTPUT=$root2.v2.3.bam.pisorted SORT_ORDER=coordinate
java -Xmx4g -jar /usr/local/bioinfo/src/picard-tools/current/SortSam.jar INPUT=$root3.v2.3.bam OUTPUT=$root3.v2.3.bam.pisorted SORT_ORDER=coordinate
java -Xmx4g -jar /usr/local/bioinfo/src/picard-tools/current/SortSam.jar INPUT=$root4.v2.3.bam OUTPUT=$root4.v2.3.bam.pisorted SORT_ORDER=coordinate

rm $root1.v2.3.bam $root2.v2.3.bam $root3.v2.3.bam $root4.v2.3.bam
mail -s "BHW AOSW mapping and sorting is done - dedup starts" "tbleroy@pierroton.inra.fr"

# picard MarkDuplicates
java -Xmx4g -jar /usr/local/bioinfo/src/picard-tools/current/MarkDuplicates.jar INPUT=$root1.v2.3.bam.pisorted OUTPUT=$root1.v2.3.bam.pisorted.dedup METRICS_FILE=$root1.v2.3.bam.pisorted.dedup.metrix.txt VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 REMOVE_DUPLICATES=true
java -Xmx4g -jar /usr/local/bioinfo/src/picard-tools/current/MarkDuplicates.jar INPUT=$root2.v2.3.bam.pisorted OUTPUT=$root2.v2.3.bam.pisorted.dedup METRICS_FILE=$root2.v2.3.bam.pisorted.dedup.metrix.txt VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 REMOVE_DUPLICATES=true
java -Xmx4g -jar /usr/local/bioinfo/src/picard-tools/current/MarkDuplicates.jar INPUT=$root3.v2.3.bam.pisorted OUTPUT=$root3.v2.3.bam.pisorted.dedup METRICS_FILE=$root3.v2.3.bam.pisorted.dedup.metrix.txt VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 REMOVE_DUPLICATES=true
java -Xmx4g -jar /usr/local/bioinfo/src/picard-tools/current/MarkDuplicates.jar INPUT=$root4.v2.3.bam.pisorted OUTPUT=$root4.v2.3.bam.pisorted.dedup METRICS_FILE=$root4.v2.3.bam.pisorted.dedup.metrix.txt VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 REMOVE_DUPLICATES=true
mail -s "BHW AOSW dedup is done - dedup starts" "tbleroy@pierroton.inra.fr"
rm $root1.v2.3.bam.pisorted $root2.v2.3.bam.pisorted $root3.v2.3.bam.pisorted $root4.v2.3.bam.pisorted

#samtools merge
samtools merge $root.v2.3.bam.pisorted.dedup.merged $root1.v2.3.bam.pisorted.dedup $root2.v2.3.bam.pisorted.dedup $root3.v2.3.bam.pisorted.dedup $root4.v2.3.bam.pisorted.dedup

