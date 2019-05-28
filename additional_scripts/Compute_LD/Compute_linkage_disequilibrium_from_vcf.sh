#!/bin/bash
#SBATCH -J "vcf"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p small
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=YOUREMAIL
#SBATCH --time=24:00:00
#SBATCH --mem=12G

input=$1 #vcffile (population specific)
pop=$2   #name of the pop for which we compute ld

mkdir pop."$pop"  #folder to store the results

#run vcftools 
vcftools --gzvcf batch.$pop.recode.vcf.gz \
       --geno-r2 \
       --max-missing 1\
       --maf 0.15 \
       --min-meanDP \
       --ld-window-bp 1000000 \
       --ld-window 1000\
       --out pop."$pop"/"$pop"

grep -v "na" pop."$pop"/"$pop".1000kb.geno.ld |grep -v "NW" >  pop."$pop"/"$pop".ld.tmp
cut -f 2,3,5 pop."$pop"/"$pop".ld.tmp  |\
        awk -f so.awk |\
        cut -f 3- >pop."$pop"/"$pop".ld_v2

rm *tmp

gzip pop."$pop"/"$pop".ld_v2
gzip pop."$pop"/"$pop".1000kb.geno.ld
