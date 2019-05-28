# compile info from the 3 vcf
python script_parser_vcf.py merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS Osativa_merged_joint_bwa_mem_mdup_raw.filtered.vcf Omeridionalis_merged_joint_bwa_mem_mdup_raw.filtered.vcf > merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withOsativa_Omeridionalis

# generate summary
echo "chr   pos    48indAll1  48indAll2  sativa_all1 sativa_all2    meridionalis_all1  meridionalis_all2" > merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withOsativa_Omeridionalis.summaryalleles
awk '{print $1"	"$2"	"$4"	"$5"	"$60"	"$61"	"$71"	"$72}' merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withOsativa_Omeridionalis >> merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withOsativa_Omeridionalis.summaryalleles
#chr    pos 48indAll1   48indAll2  sativa_all1  sativa_all2 meridionalis_all1   meridionalis_all2
#1	1248	G	A	G	.	G	.
#1	1249	A	C	A	.	A	C
#1	1266	G	A	G	.	G	.
#1	1277	T	C	T	.	T	C

# generate counts of each genotypes (double homozygous ancestral , heterozygous, double homozygous derived)
python script_ancestral_derived_counts.py merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withOsativa_Omeridionalis.summaryalleles > merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withOsativa_Omeridionalis.summaryalleles.counts

# generate derived allele frequencies
echo "chr   pos ancestral   derived missingrateOb   missingrateOg   derivedallfreqOb    derivedallfreqOg" > merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withOsativa_Omeridionalis.summaryalleles.counts.derivedallfreq
python script_compute_derivedallelefreq.py merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withOsativa_Omeridionalis.summaryalleles.counts >> merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withOsativa_Omeridionalis.summaryalleles.counts.derivedallfreq
