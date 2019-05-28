#$ -M thibault.leroy@umontpellier.fr
#$ -m abe

# download the Asian rice reference genome
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-42/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna_sm.toplevel.fa.gz
gunzip Oryza_sativa.IRGSP-1.0.dna_sm.toplevel.fa.gz
# specify the file path to the reference genome
refile=$(echo "/sandbox/users/tleroy/AfricanRice/Oryza_sativa/genome/Oryza_sativa.IRGSP-1.0.dna_sm.toplevel.fa")
# BWA
/sandbox/users/tleroy/Software/bwa-0.7.17/bwa index -a bwtsw $refile
# Samtools
/sandbox/users/tleroy/Software/samtools-1.3.1/samtools faidx $refile
# picard
module load java
java -Xmx4g  -jar /sandbox/users/tleroy/Software/picard-tools-1.140/picard.jar  CreateSequenceDictionary REFERENCE=$refile  OUTPUT=$refile.dict 
