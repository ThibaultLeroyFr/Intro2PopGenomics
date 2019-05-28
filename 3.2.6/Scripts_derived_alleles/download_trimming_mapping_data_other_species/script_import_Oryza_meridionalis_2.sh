# TL - 250219
#$ -M thibault.leroy@umontpellier.fr
#$ -m abe
# import SRA https://www.ebi.ac.uk/ena/data/view/PRJEB21312

# set pwd !
cd AfricanRice/Oryza_meridionalis

cd SRR1528445
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR152/005/SRR1528445/SRR1528445_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR152/005/SRR1528445/SRR1528445_2.fastq.gz

cd ..

cd SRR1528446
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR152/006/SRR1528446/SRR1528446_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR152/006/SRR1528446/SRR1528446_2.fastq.gz

cd ..

