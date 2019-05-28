#$ -M thibault.leroy@umontpellier.fr
#$ -m abe
#$ -q max-1m.q
#$ -pe make 4

### INFO
genus=$(echo "Oryza")
species=$(echo "glaberrima")
listacc=$(echo "ERR2009027 ERR2009041")
refile=$(echo "/sandbox/users/tleroy/AfricanRice/Oryza_sativa/genome/Oryza_sativa.IRGSP-1.0.dna_sm.toplevel.fa" )  #reference file (need to be indexed => script_index.sh) ! 

pathtodata=$(echo "/sandbox/users/tleroy/AfricanRice/"$genus"_"$species"/") # the repertory containing all individus 
pathtoscripts=$(echo "/sandbox/users/tleroy/AfricanRice/scripts/PipelineMappingCalling/") # set path for 1_mapping.sh and 2_snpindel_callingGVCF.sh 
# Please change file path in 1_mapping.sh and in 2_snpindel_callingGVCF.sh ! 

module load java # load java if needed for your cluster # GATK requires java8 !

### commandline: 

cd $pathtodata
for j in $listacc; do
    cd $j
    # CMD bash 1_mapping.sh RawData_Directory File_Trim_Paired_1 File_Trim_Paired_2 File_Trim_Unpaired_1 File_Trim_Unpaired_2 Reference_Genome Number_of_CPU_to_use
    bash $pathtoscripts/1_mapping.sh $j $pathtodata/$j/Trimmomatic *_1_cleaned.fastq.gz *_2_cleaned.fastq.gz *_1_cleaned_unpaired.fastq.gz *_2_cleaned.unpaired.fastq.gz $refile 4
    # CMD bash 2_snpindel_callingGVCF.sh IDname Reference_Genome output_directory Number_of_CPU_to_use 
    bash $pathtoscripts/2_snpindel_callingGVCF.sh $j $refile $pathtodata/$j 4
    cd ..
done
