# TL - 280319
#$ -m abe
# -q unlimitq 
#$ -l mem=150G
#$ -l h_vmem=170G

# define file & directory names (full path)
vcffile=$(echo "/bigvol/benoit/Analyses/Temp_Tibo/merged_joint_bwa_mem_mdup_raw.filtered.vcf")
outputdirscaffolds=$(echo "/bigvol/benoit/Analyses/Temp_Tibo/AfricanRice_fasta_files_scaffold")
cutoffqualitybases=$(echo "20") # minimum illumina quality at the base position >= 20
cutoffcovmin=$(echo "3") # note that this is an absolute cutoff and can be adjusted for each individual to higher values based on the distrib of coverage over the genome
cutoffcovmax=$(echo "50") # absolute threshold that can be adjusted to lower values based on distrib of coverage, for more details see Leroy et al.
outprefix=$(echo "AfricanRice_Seq")

# MAIN
### VCF2Fasta (GENERATE ALIGNED FASTA SEQUENCES 2xNbInd for each scaffold)
if [ -d "$outputdirscaffolds" ]; then
    rm $outputdirscaffolds/*.fst
else
    mkdir "$outputdirscaffolds"
fi
bash ./script_VCF2Fasta.sh $vcffile $outputdirscaffolds $cutoffqualitybases $cutoffcovmin $cutoffcovmax


