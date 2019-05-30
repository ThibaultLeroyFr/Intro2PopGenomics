# TL - 280519


# Usage: i_baypass -npop [nb_pops] -gfile [INFILE] -efile [FILE_with_evironmental] -poolsizefile [FILE_with_popsizes] -scalecov [+ parameters related to Markov chains, read the BayPass's manual]

#module load compiler/intel-2013.3.174
/usr/local/bioinfo/src/BayPass/baypass_2.1/bin/i_baypass -npop 18 -gfile BNP-BHW_18pops_v2.3_freqSNP_mincount10_filterd_biallelic_sites_counts.1.genobaypass -efile ../baypass/BNP-BHW_18pops_Toc_Pluvio.txt -poolsizefile ../baypass/BNP-BHW_18pops_popsizes.txt -scalecov -npilot 20 -thin 50 -nval 2000 -burnin 10000 -outprefix BNP-BHW_18pops_v2.3_freqSNP_mincount10_filterd_biallelic_sites_counts.1.baMeanAnTocMeanPluvio


## infiles:
# efile (populations in the same order than in the allele count matrix, line1 = temperature, line2 = precipitation sums...)
# 11.7793666667   12.0508166667   9.9597241667    9.4957683333    10.6471333333   10.59237        10.6279975      10.2194291667   8.2998316667    8.4987310833    12.32547        12.2439791667   9.2739483333    9.1595116667    7.3448365       7.1178304167    5.1583420833    6.5752054167
#786.211 791.5118        1361.6018       750.5117        698.0918        801.4064        742.6465        797.1135        635.0625        596.9492        901.0103        978.903 914.0556        933.5883        1030.5906       1016.2667       1163.8395       982.4829

# popsizes (pop in the same order than in the the allee count matrix & efile)
# 50      50      50      50      50      50      50      50      50      50      40      40      40      40      40      40      20      36



###  How to perform BayPass analyses using large sets of SNPs (> 1 million SNPs)
# Bayesian methods are highly efficient but are computationally demanding (~exponentially increasing with the number of SNPs)
# To decrease the computation time, the best thing to do is to generate reduced datasets (e.g. < 1 million SNPs for BayPass). 
 
# To do so, two options: 

# Option 1/ if you use the synchronized file (popoolation2), it is easier to use the R package poolfstat (function) to subsample the data.
# (see generate_inputfile.R for details)
#
# then, perform the computations under BayPass


# Option 2/ if your data are available on another format, the input file is a simple matrix of allele counts:
# counts_all1_pop1	counts_all2_pop1	counts_all1_pop2	counts_all2_pop2	...
#106     66      101     70    ...
#
# Based on this matrix, it is possible to get one line every 100 line (or 10 or 1000 depending of the total number of SNPs), by doing:
#for i in {1..100}; do
#        sed -n "$i~100p"  [infile] > [infile].$i
#done
#
# Then, it is possible to perform the BayPass computation for each subset.
#
#


# In both cases, final BayPass output files need to be merged (it is only possible if the inferred pop structure (omega matrix) is roughly nearly identically inferred for the different subsets... )
# paste require to be sure that the number of lines in the two sets are identical (position files = snpdet for poolfstat)
## compil positions + XtX
#for i in {1..100}; do filextx=$(echo "BNP-BHW_18pops_v2.3_freqSNP_mincount10_filterd_biallelic_sites_counts.""$i"".LUinCG_unscaled_summary_pi_xtx.out"); file=$(echo "BNP-BHW_18pops_v2.3_freqSNP_mincount10_filterd_biallelic_sites_counts.""$i"".positions"); echo "### $i ###"; grep -v "MRK" $filextx > bidon; paste $file bidon >> BNP-BHW_18pops_v2.3_freqSNP_mincount10_filterd_biallelic_sites_counts.LUinCG_unscaled_summary_pi_xtx.out.merged.withpositions; done

# Compil BF
#for i in {1..100}; do fileBF=$(echo "BNP-BHW_18pops_v2.3_freqSNP_mincount10_filterd_biallelic_sites_counts.""$i"".LUinCG_unscaled_summary_betai_reg.out");  echo "### $i ###"; awk '$1 == 1 {print $0}' $fileBF >> BNP-BHW_18pops_v2.3_freqSNP_mincount10_filterd_biallelic_sites_counts.merged.LUinCG_unscaled_summary_betai_reg.out.COVAR1; done
