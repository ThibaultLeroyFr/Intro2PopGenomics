# TL - 280319
#/bigvol/benoit/Analyses/Mapping/F_albicollis/joint_albicollis_bwa_mem_mdup_filtered.vcf
file=$1
outputdir=$2
cutoffquality=$3 #20
cutoffcovmin=$4
cutoffcovmax=$5

infiledir=$(dirname "${file}")
infilename=$(basename "${file}")

# create the outputfile, report an error if this directory already exists
if [ ! -d $outputdir  ]; then
    mkdir -p $outputdir
else
#    echo "[ERROR] the followig output directory $outputdir already exists, premature ending of the script"
#    exit
    continue
fi

# check the path, unzip the vcf file if found or report error
if [ -f $file ]; then
    echo "[INFO] starting initial quality check"
elif [ -f $file.zip ]; then #if a zip file of the vcf is found, unzip it
    echo "[INFO] compressed file found, starting of unzipping"
    unzip -d $infiledir $file.zip
    echo "[INFO] end of the unzipping step"
    echo "[INFO] starting intial quality check"
else
    echo "[ERROR] File $infilename at $infiledir not found, premature ending of the script"
    exit
fi

# extract quality at both SNP and individuel level
#grep -v "#" $file | grep -v "LowQual" | awk '$7 != "." {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$8"\t"$7}' > $file.qual
#python ./scripts/FiltrationDepthCoverageQualityJointVCF_PASSLowQ_variantonly.py $file.qual > $file.qualsummary
#python ./scripts/computeDepthCoverageQualityJointVCF_variantonly.py $file > $file.qualsummary.ind
#awk '$4 >= 0 {print $0}' $file.qualsummary.ind > $file.qualsummary.ind.clean #rm non-numeric values if exists

# compute 99CI confidence intervals for the coverage at the SNP level and at the individual level
#sed "s;~/Collab/Montpellier/Colivacea_CovQual_joinVcf/Colivacea_FiltrationQualityjointVCF_variantonly_allsites_060218.PROV.30pc.sed.withoutNA;$file.qualsummary;g" ./scripts/script_coverage_R_gvcf.R | sed "s;C_olivacea_summary_calling_newFilteredGVCF_060218.pdf;$file.qualsummary.pdf;g" | sed "s;summary_statistics_coverage;$file.99CI;g" | sed "s;~/Collab/Montpellier/Colivacea_CovQual_joinVcf/Colivacea_bidon;$file.qualsummary.ind.clean;g" >  ./scripts/script_coverage_R_gvcf.$infilename.tmp
#echo "[INFO] starting R computations"
#Rscript ./scripts/script_coverage_R_gvcf.$infilename.tmp
 #echo "[INFO] ending R computations"
#head $file.99CI.txt
#head $file.99CI_ind.txt
#rm script_coverage_R_gvcf.R.tmp

### generating fasta files
echo "[INFO] starting to create fasta files from VCFs"
cd $outputdir
python ../VCF2Fasta_fast.py -q $cutoffquality -m $cutoffcovmin -M $cutoffcovmax -f $file > $outputdir/Outputs_VCF2Fasta.txt
echo "[INFO] fasta files created"
time=$(date)
echo "[INFO] Computations performed succesfully - this script finished at $time"
echo "############################################################################################"





