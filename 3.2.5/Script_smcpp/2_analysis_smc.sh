#!/bin/bash
#SBATCH -J "smcpp_cv_8e-09_69k"
#SBATCH -o log_%j
#SBATCH -c 5
#SBATCH -p small
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=YOUREMAIL
#SBATCH --time=22:30:00
#SBATCH --mem=06G

# Move to directory where job was submitted
#cd $SLURM_SUBMIT_DIR

TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

pop_arg=$1
if [[ -z "$pop_arg" ]]
then
    echo "error please provide pop_name"
    exit
fi

pop=$(echo "$pop_arg")

outfolder="10.analysis_estimate"
mkdir -p "$outfolder"/"$pop"

#thin="--thinning 1609" #se value for thinning
#t1="--t1 50" #set the min number of generation
#tk="--tK 250000" #set the max nb generation
#ft="--ftol 0.01"
#core="--cores 5"
#out="-o 03.analysis_v1"
#rp="--regularization-penalty 5"
#pe="--polarization-error 0.5"

#smc++ estimate "$core" "$out" "$thin" "$t1" "$tk" "$ft" "$rp" "$pe" 0.5e-8 02.out/"$pop"/chr.*.smc.gz
smc++ cv --cores 5 \
	--out "$outfolder"/$pop \
	#--thinning 1609 \
	#--timepoints 50 500000 \
	----Nmax 1000000\
	#--ftol 0.0001 \
	#--regularization-penalty 6\
	#--polarization-error 0.5 \
	--knots 25\
	6.5e-9 02.out/"$pop"/chr.*.smc.gz #'2>&1' '>>' 10-log/"$TIMESTAMP"_anlaysis"$pop".log
