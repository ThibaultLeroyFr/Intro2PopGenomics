# TL - 280519
# Rscript generate_inputfile.R

library("poolfstat") # to install: install.packages("poolfstat")

setwd("/home/thibaultleroy/Collab/Sessile_Reanalyse")
# read the synchronized file (popoolation2)
# details about the data (poolsizes = number of haplotypes = 2*nb_diploid_individuals; populations names; minimum number of reads supporting the alternate alleles, minimum and maximum coverage for the pop, minor allele frequency etc...
pooldata=popsync2pooldata(sync.file="/media/thibaultleroy/Seagate Expansion Drive/projet_BNP-BHW/PileupSync/BNP-BHW_18pops_v2.3.pileup.sync.1",poolsizes=c(50,50,50,50,50,44,50,50,50,50,40,40,40,40,40,40,20,36),poolnames=c("9","97","124","204","217","218","219","233","253","256","L1","O1","L8","O8","O12","L12","O16","L16"),min.rc=10,min.cov.per.pool = 50, max.cov.per.pool = 250, min.maf=0.02,noindel=TRUE)

# Then convert the data
# writing.dir=getwd() = print the files in the current directory
# subsampling size = create two sets of half SNPs (subsamplesize=1 => only convert, no subsampling)
# subsampling is very important to reduce the computation time. To my opinion, the number of SNPs need to be 50,000 - 1,000,000 for each subset (too few loci = inferred differences in omega matrix between BayPass runs; too many loci = ) 
pooldata2genobaypass(pooldata,writing.dir=getwd(),prefix="",subsamplesize=200000,subsamplingmethod="thinning")

# Expected output files
# "genobaypass.sub" output files (allele counts, for each pop counts_allele1 counts_allele2 ... ) => input for baypass
#101 5 138 5 110 14 123 4 146 5 140 9 133 8 122 8 174 11 178 7 135 5 103 3 128 5 101 2 124 10 120 5 133 6 123 9
#65 11 72 24 40 18 67 20 90 6 78 17 72 12 62 9 141 14 118 24 73 8 70 8 72 19 66 12 52 21 43 20 59 20 58 15

# "snpdet.sub" output files corresponding to SNP info
#Sc0000000 8300 c T C
#Sc0000000 72750 G T G

