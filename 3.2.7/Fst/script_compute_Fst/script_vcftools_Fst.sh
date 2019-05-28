# TL - 280519
# this script use vcftools to compute Weir & Cockerham Fst values
# usage: vcftools --vcf [file] --weir-fst-pop [file_ID_ind_pop1] --weir-fst-pop [file_ID_ind_pop2] --fst-window-size [window_size] --fst-window-step [step_btw_2_computations]

# To compute Fst on sliding windows (add as many "--weir-fst-pop" arguments as required, in case you work on more than 2 populations):
vcftools --vcf merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withheader --weir-fst-pop Oryza_barthii_individuals --weir-fst-pop Oryza_glaberrima_individuals --fst-window-size 100000 --fst-window-step 100000


# Expected output:
#CHROM   BIN_START       BIN_END N_VARIANTS      WEIGHTED_FST    MEAN_FST
#1       1       100000  884     0.0705145       0.0304594
#1       100001  200000  839     0.0526089       0.0226013
#1       200001  300000  1025    0.150095        0.0666252
#...


# To compute Fst on a per-SNP basis (--fst-window-size 1 --fst-window-step 1 by default, add as many "--weir-fst-pop" arguments as required, in case you work on more than 2 populations):
vcftools --vcf merged_joint_bwa_mem_mdup_raw.filtered.vcf.PASS.withheader --weir-fst-pop Oryza_barthii_individuals --weir-fst-pop Oryza_glaberrima_individuals

#Expected output:
#CHROM   POS     WEIR_AND_COCKERHAM_FST
#1       1248    -nan
#1       1249    -nan
#1       1266    -nan
#1       1277    -nan
#1       1280    -nan
#1       1299    -nan
#1       1301    -nan
#1       1352    0.00589264
# ...

# useful commands:
# sed 's/-nan/NA/g' to edit -nan to NA 
# grep -v "\-nan" to remove all lines without computed values.
