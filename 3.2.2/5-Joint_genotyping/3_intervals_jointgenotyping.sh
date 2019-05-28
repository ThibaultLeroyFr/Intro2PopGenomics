assembly=$1  # genome assembly
intervals=$2   # expected nb intervals (number of intervals created can be slightly lower or higher)

module load bioinfo/Java8 # very important for the GATK joint calling

# tmp dir
assemblyname=$(basename $assembly)
mkdir tmp_vcf_$assemblyname

# create intervals
python /home/tleroy/work2/Collab/AfricanRice/script_scaff_length.py $assembly > ./tmp_vcf_$assemblyname/$assemblyname.scaffsize
python /home/tleroy/work2/Collab/AfricanRice/createintervalsfromscaffsize.py ./tmp_vcf_$assemblyname/$assemblyname.scaffsize $intervals
mv scatter*.intervals ./tmp_vcf_$assemblyname/
cd ./tmp_vcf_$assemblyname/

# loop over intervals
for i in scatter*.intervals; do
	#echo "time java -jar /media/bigvol/benoit/software/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -T GenotypeGVCFs -nt 1 -R $assembly -o $i.joint_bwa_mem_mdup_raw.vcf -V ../gcfFileList.list  -L $i" > lanceur_interval.sh
	echo "#$ -m abe" > lanceur_interval.$i.sh
	echo "#$ -l mem=20G" >> lanceur_interval.$i.sh
        echo "#$ -l h_vmem=22G" >> lanceur_interval.$i.sh
    	echo "module load bioinfo/Java8" >> lanceur_interval.$i.sh
	echo "java -Xmx4g -jar ~/work/Software/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -T GenotypeGVCFs -nt 1 -R $assembly  -o $i.joint_bwa_mem_mdup_raw.vcf -V ../gcfFileList.list  -L $i --includeNonVariantSites" >> lanceur_interval.$i.sh
	qsub lanceur_interval.$i.sh # qsub or sbatch for genotoul
    sleep 1
done
# merging at the end
#for i in {1..$intervals}; do file=$(echo "scatter""$i"".intervals.joint_bwa_mem_mdup_raw.vcf"); if [ $i == 1 ]; then cp $file merged_joint_bwa_mem_mdup_raw.vcf; else grep -v "#" $file >> merged_joint_bwa_mem_mdup_raw.vcf; fi; done
