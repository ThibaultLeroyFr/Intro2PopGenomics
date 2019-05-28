#!/bin/bash
#SBATCH -J "vcfToSmcc"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p small
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=YOUREMAIL
#SBATCH --time=22:30:00
#SBATCH --mem=07G

# Move to directory where job was submitted
#cd $SLURM_SUBMIT_DIR
source ~/software/environments/my_env/bin/activate
#activate python3 environment:
source activate python3 #this has been seted before

#load variables
input="output_PASS_only"
input_folder="01.data"

pop_arg=$1
if [[ -z "$pop_arg" ]]
then
       echo "error please provide pop_name"
       exit
fi

pop=$(echo "$pop_arg" |sed 's/pop.//g')
ind=$(less "$pop_arg" |perl -pe 's/\n/,/g' |sed 's/\,$//g')

mkdir 02.out/"$pop" 2>/dev/null

for i in $(cat list_NC_2 ) ;
do
    smc++ vcf2smc --cores 5\
    "${input_folder}"/"${input}"."$i".recode.vcf.gz \
    02.out/"$pop"/chr."$i".smc.gz \
    "$i" \
    "$pop":"$ind"
done


