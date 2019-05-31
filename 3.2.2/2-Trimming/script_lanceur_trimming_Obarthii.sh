# TL - 050319
#$ -m abe
# -pe parallel_environment 10
module load java
for j in ERR2008849 ERR2008850 ERR2008851 ERR2008852 ERR2008853 ERR2008855 ERR2008856 ERR2008857 ERR2008858 ERR2008859 ERR2008860 ERR2008861 ERR2008862 ERR2008863 ERR2008884 ERR2008894 ERR2008910 ERR2008914 ERR2009084 ERR2009085 ERR2009086 ERR2009090 ERR2009096; do
        cd /sandbox/users/tleroy/AfricanRice/Oryza_barthii/$j
        mkdir Trimmomatic
        seqlines=$(ls | grep -v "MD5" | grep -v "Trimmomatic" | grep -v "total" | sed 's/_1.fastq.gz//g' | sed 's/_2.fastq.gz//g' | sort | uniq)
        for i in $seqlines; do
            acc=$(echo "/sandbox/users/tleroy/AfricanRice/Oryza_barthii/ERR2008849/$j/$i")
            outacc=$(echo "/sandbox/users/tleroy/AfricanRice/Oryza_barthii/ERR2008849/$j/Trimmomatic/$i")
            java -Xmx4g -jar /sandbox/users/tleroy/Software/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 1 -phred33 "$acc"_1.fastq.gz "$acc"_2.fastq.gz "$outacc"_1_cleaned.fastq.gz "$outacc"_1_cleaned_unpaired.fastq.gz "$outacc"_2_cleaned.fastq.gz "$outacc"_2_cleaned_unpaired.fastq.gz ILLUMINACLIP:/sandbox/users/tleroy/Software/Trimmomatic-0.33/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
    done
    #rm "$acc"_1.fastq.gz "$acc"_2.fastq.gz
done
