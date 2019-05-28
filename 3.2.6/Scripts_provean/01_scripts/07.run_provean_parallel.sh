#!/bin/bash
#SBATCH --account=rrg-blouis
#SBATCH --time=48:00:00
#SBATCH --job-name=abc
#SBATCH --output=abc-%J.out
##SBATCH --array=1-33%33
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
##SBATCH --gres=cpu:32
##SBATCH --mail-user=yourmail
##PBS -l nodes=1:ppn=8
##SBATCH --mail-type=EA 

# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

# External args
list_NC=$1
if [ -z "$list_NC" ]
then
	echo "error ! Need chromosome folder"
	echo exit
fi

if [ ! -d "08_provean_scores" ]
then
	echo "creating folder"
	mkdir 08_provean_scores"
fi
if [ ! -d "07_provean_match" ]
then
	echo "creating folder"
	mkdir 07_provean_match"
fi

mkdir -p 08_provean_scores/"${list_NC}"
mkdir -p 07_provean_match/"${list_NC}"

# Global variables
input="$list_NC"

# Launching provean in parallel
find /home/quentin/scratch/24.rice/06_provean/"$input" | \
	grep "\.fasta$" | \
	parallel -j 32 /home/quentin/scratch/24.rice/01_scripts/08_run_provean_iteration.sh {} "$input"
