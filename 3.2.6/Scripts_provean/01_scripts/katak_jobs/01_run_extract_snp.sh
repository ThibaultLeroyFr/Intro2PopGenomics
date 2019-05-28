#!/bin/bash
#SBATCH -J "seq_extract"
#SBATCH -o log_%j
#SBATCH -c 1
#SBATCH -p ibismax
#SBATCH -A ibismax
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=YOUREMAIL
#SBATCH --time=4-30:00
#SBATCH --mem=08G

# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

# Load cutadapt module
FASTA=$1
WANTED_LOCI=$2
LENGTH=$3 

if [[ -z $FASTA ]]
then
    echo "error need fasta input file"
    exit
fi

if [[ -z $WANTED_LOCI ]]
then
    echo "error need loci input file"
    exit
fi

if [[ -z $FASTA ]]
then
    echo "error need length "
    exit
fi

# Reduce number of CPUs here and above if you have less than 4 chips/lanes
./01_scripts/01_extract_snp_variants_with_flanking.py $FASTA $WANTED_LOCI $LENGTH ${WANTED_LOCI%}.fasta
