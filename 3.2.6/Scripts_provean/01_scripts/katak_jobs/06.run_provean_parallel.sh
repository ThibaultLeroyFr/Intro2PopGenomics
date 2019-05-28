#!/bin/bash
#SBATCH -J "wgstest"
#SBATCH -o log_%j
#SBATCH -c 10
##SBATCH -p low-suspend
#SBATCH -p ibismax
#SBATCH -A ibismax
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=YOUREMAIL
#SBATCH --time=20-00:00
#SBATCH --mem=12G

# Move to directory where job was submitted
#cd $SLURM_SUBMIT_DIR

# External args
list_NC=$1
if [ -z "$list_NC" ]
then
	echo "error ! Need chromosome folder"
	echo exit
fi

# Global variables
input="$list_NC"

# Launching provean in parallel
/bin/bash: q : commande introuvable
