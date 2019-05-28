# simple faststructure analyses for K=1 to K=10
for i in {1..10}; do
	echo "starting simulation for K=$i"
	python structure.py -K $i --input merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withoutNA.3pc.withheader.bed --output merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withoutNA.3pc.withheader.K --full --cv 10 --format bed
done
