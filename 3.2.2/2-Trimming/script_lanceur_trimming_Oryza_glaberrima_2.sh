# TL - 050319
#$ -M thibault.leroy@umontpellier.fr
#$ -m abe
# -pe parallel_environment 10
module load java
for j in ERR2008923 ERR2008927 ERR2008934; do
        cd /sandbox/users/tleroy/AfricanRice/Oryza_glaberrima/$j
        mkdir Trimmomatic
        seqlines=$(ls | grep -v "MD5" | grep -v "Trimmomatic" | grep -v "total" | sed 's/_1.fastq.gz//g' | sed 's/_2.fastq.gz//g' | sort | uniq)
        for i in $seqlines; do
            acc=$(echo "/sandbox/users/tleroy/AfricanRice/Oryza_glaberrima/$j/$i")
            outacc=$(echo "/sandbox/users/tleroy/AfricanRice/Oryza_glaberrima/$j/Trimmomatic/$i")
            java -Xmx4g -jar /sandbox/users/tleroy/Software/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 1 -phred33 "$acc"_1.fastq.gz "$acc"_2.fastq.gz "$outacc"_1_cleaned.fastq.gz "$outacc"_1_cleaned_unpaired.fastq.gz "$outacc"_2_cleaned.fastq.gz "$outacc"_2_cleaned_unpaired.fastq.gz ILLUMINACLIP:/sandbox/users/tleroy/Software/Trimmomatic-0.33/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
    done
    rm "$acc"_1.fastq.gz "$acc"_2.fastq.gz
done


