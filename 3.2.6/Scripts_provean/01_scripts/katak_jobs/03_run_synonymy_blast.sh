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
INPUT_FASTA=$1
DATABASE=$2

if [[ -z $INPUT_FASTA ]]
then
    echo "error need fasta input file"
    exit
fi

if [[ -z $DATABASE ]]
then
    echo "error need blast database from protein seq"
    exit
fi

./01_scripts/03_synonymy_blast.sh  "$INPUT_FASTA" "$DATABASE"
