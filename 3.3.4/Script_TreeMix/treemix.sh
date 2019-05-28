# TL - 280519

# IMPORTANT: the input file is a simple allele count matrix and need to be gzipped 
# Assuming that the input file is: "bidon.infiletreemix.gz", [infile] = "bidon.infiletreemix" (without the ".gz" in the file name)

# Expected infile (i.e. zmore bidon.infiletreemix.gz ) 
#Line1: pop_name_1	pop_name_2	popname_3 ... *n cols (each column is a population)
#Line2 (SNP#1): count_all1pop1,count_all2pop1	count_all1pop2,count_all2pop2	count_all1pop2,count_all2pop2...
# ...
# n rows (each row is a SNP)


# Treemix only requires single-command lines.
# Two key parameters:
# k (grouping of neighboring SNPs to take into account LD) = 1000
# m (number of migration nodes assumed) e.g. m=0 (default value)

# e.g. m=0
./treemix-1.13/src/treemix -i [infile] -k 1000 -o [outfile].m.0
# m1
./treemix-1.13/src/treemix -i [infile] -k 1000 -m 1 -o [outfile].m.1
# ...


# or using a simple loop:
for i in {1..10}; do 
	./treemix-1.13/src/treemix -i [infile] -k 1000 -m $i -o [outfile].m.$i
done

