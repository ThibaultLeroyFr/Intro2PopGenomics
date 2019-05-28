#!/bin/bash
##SBATCH -J "smcpp"
##SBATCH -o log_%j
##SBATCH -c 5
##SBATCH -p small
##SBATCH --mail-type=FAIL
##SBATCH --mail-user=YOUREMAIL
##SBATCH --time=22:30:00
##SBATCH --mem=08G

# Move to directory where job was submitted
#cd $SLURM_SUBMIT_DIR

#TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
#echo $arg
infile=$1
mkdir "$infile"_plot
smc++ plot "$infile"_plot/plot_generation.pdf -g 1 -c "$infile"/*/model.final.json  
smc++ plot "$infile"_plot/plot_coalescent.pdf --csv  "$infile"/*/model.final.json   
exit
smc++ plot "$infile"_plot/plot_generation_knots.pdf -g 4 --knots  "$infile"/*/model.final.json   
exit
smc++ plot 04.plot/plot_generation.pdf "$arg" -g 4 -c 
smc++ plot 04.plot/plot_coalescent.pdf "$arg" -g 4 
smc++ plot 04.plot/plot_generation_knots.pdf "$arg" -g 4  
