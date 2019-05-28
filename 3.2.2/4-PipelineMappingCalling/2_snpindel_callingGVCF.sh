#! /usr/bin/env bash
## script: 'snpindel_callingGVCF_raw.sh'

## your GATK sample
f=$1 # IDname
refgen=$2 # reference genome
infolder=$3 #infolder same than outfolder
nthr=$4 # nb of cpu to use
echo $3

tmpdir=$(echo "/tmp")

# Requires java8 !
java -Xmx4g -Djava.io.tmpdir=$tmpdir -jar /sandbox/users/tleroy/Software/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -nct ${nthr} -T HaplotypeCaller \
   -R ${refgen} -I ${infolder}/${f}_bwa_mem_mdup.bam \
   --variant_index_type LINEAR \
   --heterozygosity  0.002 \
   --variant_index_parameter 128000 \
   --emitRefConfidence GVCF  \
   --genotyping_mode DISCOVERY  \
   -o ${infolder}/${f}_bwa_mem_mdup_raw.g.vcf
   
   
echo "# finished Raw Variant Calling GVCF for sample ${f} with GATK"

