#$ -M thibault.leroy@umontpellier.fr
#$ -m abe

# download the Q. robur genome (a closely-related European white oak species)
https://urgi.versailles.inra.fr/download/oak/Qrob_PM1N.fa.gz
gunzip Qrob_PM1N.fa.gz
# specify the file path to the reference genome
refile=$(echo "/sandbox/users/tleroy/SessileOak/Qrob_PM1N.fa")
# BWA
/sandbox/users/tleroy/Software/bwa-0.7.17/bwa index -a bwtsw $refile
# Bowtie2 index
#Usage: bowtie2-build [options]* <reference_in> <bt2_index_base>
bowtie2-build Qrob_PM1N.fa Qrob_PM1N
# Samtools
/sandbox/users/tleroy/Software/samtools-1.3.1/samtools faidx $refile
# picard
module load java
java -Xmx4g  -jar /sandbox/users/tleroy/Software/picard-tools-1.140/picard.jar  CreateSequenceDictionary REFERENCE=$refile  OUTPUT=$refile.dict 

