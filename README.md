# Introduction to population genomics
This github directory hosts all scripts used to perform the analyses shown in in the book chapter of Leroy & Rougemont (2019).

>*Leroy, T & Rougemont, Q. 2019. Population genetics analysis methods: from population structure inferences to genome scans for selection. In: Molecular Plant Taxonomy (Springer)*

Please send an email to <a href="mailto:thibault.leroy@univie.ac.at;quentinrougemont@orange.fr?subject=[Intro2popGenomics-Github]">us </a> for questions regarding these scripts or for full-text requests.

## Foreword

### Table of contents:<br/>
*The github repository follow the same organisation as in the book chapter.*<br/>
**Individual data (African rice analyses, 3.2)**<br/>
*3.2.2 : From raw data to VCF*<br/>
*3.2.3 - Population structure*<br/>
*3.2.4 - Nucleotide diversity*<br/>
*3.2.5 - Historical population sizes*<br/>
*3.2.6 - Deleterious mutation load*<br/>
*3.2.7 - Fst & genome scans*<br/>
**Pool-seq data (Sessile oak analyses, 3.3)** <br/>
*3.3.1 - Pool-seq vs. individual data* <br/>
*3.3.3 - From raw data to allele counts* <br/>
*3.3.4 - Population splits & mixtures* <br/>
*3.3.5 - Fst estimates* <br/>
*3.3.6 - Genome Scan for Selection* <br/>
*3.3.7 - Genotype-Environment associations* <br/>


### Important note
Our scripts are not standalone executables. Quite the contrary, these scripts (deliberately) require some simple edits to adjust to your data. Editing script is probably the best way to learn how a script works, to detect and correct the errors and, more broadly, to start learning how to code. So please take some time to read the scripts and, ideally, to briefly look at the software user manuals. Change the paths to files and programs to adjust the scripts to your computer architecture. 


## Overview
### Individual data - Application to the African Rice data </br>
*3.2.2 : From raw data to VCF*</br>

<pre><code>
<strong>1/Import sequencing data (./3.2.2/1-Import_RawData/)</strong>
Softwares needed: <a href="https://www.gnu.org/software/wget/">wget</a> (ftp-transfert)
<em><a href="https://www.ebi.ac.uk/ena/data/view/PRJEB21312">Here a list of all ftp addresses</a>
e.g.wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR200/000/ERR2008850/ERR2008850_1.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR200/000/ERR2008850/ERR2008850_2.fastq.gz</em>

<strong>2/Read trimming (./3.2.2/2-Trimming/)</strong>
Softwares needed: <a href="https://github.com/timflutre/trimmomatic">Trimmomatic</a>
<em>java -Xmx4g -jar ./trimmomatic-0.33.jar PE -threads 1 -phred33 "[file]_1.fastq.gz [file]_2.fastq.gz [file]_1_cleaned.fastq.gz [file]_1_cleaned_unpaired.fastq.gz [file]_2_cleaned.fastq.gz [file]_2_cleaned_unpaired.fastq.gz ILLUMINACLIP:./adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50</em>

<strong>3/Downloading & indexing a reference genome (./3.2.2/3-Download_Index_references/)</strong>
Softwares needed: <a href="https://www.gnu.org/software/wget/">wget</a> (ftp-transfert), <a href="https://sourceforge.net/projects/bio-bwa/files/">BWA</a>, <a href="http://samtools.sourceforge.net/">Samtools</a> & <a href="https://broadinstitute.github.io/picard/">Picard</a> (indexing) 
<em>Asian rice genome: wget ftp://ftp.ensemblgenomes.org/pub/plants/release-42/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna_sm.toplevel.fa.gz
BWA: bwa index -a bwtsw [fasta]
Samtools: samtools faidx [fasta]
Picard: java -Xmx4g  -jar picard.jar  CreateSequenceDictionary REFERENCE=[fasta]  OUTPUT=[fasta].dict </em>

<strong>4/Mapping & Individual Calling (./3.2.2/4-PipelineMappingCalling/)</strong>
Softwares needed: <a href="https://sourceforge.net/projects/bio-bwa/files/">BWA</a> (mapping), <a href="http://samtools.sourceforge.net/">Samtools</a> (filtering), <a href="https://broadinstitute.github.io/picard/">Picard</a> (removing duplicates) & <a href="https://software.broadinstitute.org/gatk/download/">GATK</a> (creating gVCF files)
<em>bash 1_mapping.sh [RawData_Directory] [File_Trim_Paired_1] [File_Trim_Paired_2] [File_Trim_Unpaired_1] [File_Trim_Unpaired_2] [Reference_Genome] [Number_of_CPU_to_use]
bash 2_snpindel_callingGVCF.sh [ID] [Reference_Genome] [output_directory] [Number_of_CPU_to_use] </em>

<strong>5/Joint Genotyping (./3.2.2/5-Joint_genotyping/)</strong>
<em> bash 3_intervals_jointgenotyping.sh [reference_genome] [nb_cpus]  </em>

<strong>6/Filtering variants (./3.2.2/6-Variant_Filtration/)</strong>
<em> ./VariantFiltrationVCF.py -q 2.0 -s 60.0 -m 40.0 -n -2.0 -r -2.0 -w 45000 -f [VCF] > [filtered_VCF] </em>
</pre></code>

*3.2.3 - Population structure*</br>
<pre><code>
<strong>Perform Principal Component Analysis (./3.2.3/PCA/)</strong>
R packages needed: <a href="https://bioconductor.org/packages/release/bioc/html/SNPRelate.html/">SNPRelate</a> (& <a href="https://bioconductor.org/packages/release/bioc/html/gdsfmt.html/">gdsfmt</a>)
<em>Details shown in script_PCA_from_vcf.R (Import & convert the vcf file, compute PCA & generate plots) </em>
</pre></code>

<pre><code>
<strong>Infer ancestry proportions (./3.2.3/fastStructure/)</strong>
Software needed: <a href="https://www.cog-genomics.org/plink2/">Plink</a> & <a href="https://rajanil.github.io/fastStructure/">fastStructure</a>
Generate an input file using Plink:
<em> ./plink --allow-extra-chr -allow-extra-chr --make-bed --noweb --out [VCF].bed --vcf [VCF] </em>
FastStructure inferences:
<em>./structure.py -K [nb_of_clusters] --input [VCF].bed --output [VCF].bed.K[nb_of_clusters] --full --cv [cross-validation steps] --format bed
./chooseK.py --input=[VCF].bed.K </em>
</pre></code>

*3.2.4 - Nucleotide diversity*<br/>
<pre><code>
<strong>Generate Fasta files from VCF (./3.2.4/Scripts_generate_fasta_sequences_from_vcf/)</strong>
Scripts were developed by Leroy et al. 2019 bioRxiv 505610, ver. 4 peer-reviewed and recommended by PCI Evolutionary Biology (paper: <a href="https://www.biorxiv.org/content/biorxiv/early/2019/05/24/505610.full.pdf">here</a>, scripts: <a href="https://figshare.com/s/122efbec2e3632188674">here</a>)
<em>python ../VCF2Fasta_fast.py -q [threshold_base_quality] -m [minimum_coverage_position] -M [maximum_coverage_position] -f [VCF] > [outputdir]/Outputs_VCF2Fasta.txt</em>

<strong>Compute theta, pi & Tajima's D (./3.2.4/Scripts_compute_pi_D/)</strong>
Software needed: <a href="https://figshare.com/s/122efbec2e3632188674#/articles/7484705">seq_stat</a>
<em>bash script_compute_pi_slidwin.sh [referencegenome].scafflength [VCF] [ouputdir] [size_of_sliding_window] [output_prefix]</em>
(scafflength = file containing the length of each scaffold as computed by ./3.2.2/5-Joint_genotyping/script_scaff_length.py):  

<strong>Generate a circlize plot (./3.2.4/Rscript_plot_Circlize/) </strong>
R package needed: <a href="https://cran.r-project.org/web/packages/circlize/index.html">circlize</a> (more details: <a href="https://jokergoo.github.io/circlize_book/book/">here</a>) 
<em>Details shown in script_circlize_GC_pi_ROD_Africanrice.R (Import the file containing the summary statistics & generate a circlize plot)
(the file containing the summary statistics was made available at ./3.2.4/Results_pi/Genomic_pi_ROD_010419.withoutNA.txt)</em>
</pre></code>

*3.2.5 - Historical population sizes*<br/>
Software needed: <a href="https://github.com/popgenmethods/smcpp">smc++</a>
<pre><code>
<strong>Convert the vcf to the smc++ input format (./3.2.5/1_vcf2smcpp.sh)</strong>
<em>smc++ vcf2smc --cores [nb_cpu] [input_vcf_file] [output_smc_data_files] [chr] [pop1:Ind1,Ind2,Ind3..]</em>

<strong>Perform the inference (./3.2.5/2_analysis_smc.sh)</strong>
<em>smc++ cv --cores [nb_cpu] --out [out] --Nmax [Ne_max] --knots [spline_knots_for_smoothing] [mutation_rate] [smc_data_files]</em>

<strong>Generate a plot (./3.2.5/3_smcpp_plot.sh)</strong>
<em>smc++ plot [outfile.pdf] -g 1 -c [infile_model.final.json]
(-c produces a CSV-formatted table: this file is also available ./3.2.5/Rscript_plot/plot_generation.csv)</em>
It is also possible to use the newly generated .csv file to generate the plot under the R environment, e.g. using the R package <a href="https://cran.r-project.org/web/packages/ggplot2/index.html">ggplot2</a> : <em>see ./3.2.5/Rscript_plot/script_generateplot_smcpp_220519.R for details</em>
</pre></code>


*3.2.6 - Deleterious mutation load*<br/>
<pre><code>
<strong>Counts the number of derived alleles</strong>

1/Generate VCF for the outgroup species (./3.2.6/Scripts_derived_alleles/download_trimming_mapping_data_other_species)
<em>Download the raw data for the outgroup species -> generate joint VCF / outgroup species
Same steps than in the section "3.2.2 : From raw data to VCF"</em>
    
2/Detect the ancestral state & compute allele frequencies (./3.2.6/Scripts_derived_alleles/script_vcfoutgroup_to_ancestral_derived.sh) 

In a nutshell:
Parse the 3 joint vcf (the focal species & the 2 newly obtained vcf corresponding to the 2 outgroup species) 
<em>./script_parser_vcf.py [VCF_focal_species_ONLY_PASS_variants] [VCF outgroup1] [VCF outgroup2]> [Merged_VCF_file]</em>

Then use awk '{print $X"    "Y...}' (where X and Y correspond to the columns in the [Merged_VCF_file]) to parse the data to obtain the following file format:
<em>chr    pos focal_All1   focal_All2  outgroup1_all1  outgroup1_all2 outgroup2_all1   outgroup2_all2
1	1248	G	A	G	.	G	.</em>

Detect the ancestral allele & compute allele frequencies
<em>./script_ancestral_derived_counts.py [this_infile] > [file_with_derived_counts]
./script_compute_derivedallelefreq.py [file_with_derived_counts] > [files_with_derivedallfreq]</em>
</pre></code>

<pre><code>
<strong>Deleterious variant prediction (./3.2.6/Scripts_provean/) </strong>
Softwares needed: <a href="https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download">blast+</a> & <a href="http://provean.jcvi.org/downloads.php">provean</a>
Follow each step in the order indicated (from script 01 to script 08)
</pre></code>

3.2.7 - Fst & genome scans
<pre><code>
<strong>Compute Fst from vcf (./3.2.7/script_compute_Fst/script_vcftools_Fst.sh) </strong>
Softwares needed: <a href="http://vcftools.sourceforge.net/">vcftools</a>  
<em>./vcftools --vcf [VCF] --weir-fst-pop [IDs_sp1] --weir-fst-pop [IDs_sp2] --fst-window-size [window_size_in_bp] --fst-window-step [size_window_step_in_bp]
where [IDs_sp1], [IDs_sp2] , ... correspond to a file containing a list of all individuals for sp1, sp2, ... (require the same ID as in the vcf)</em>

To compute Fst on non-overlapping sliding windows (preferred), use the same value for both [window_size_in_bp] and [size_window_step_in_bp]
To compute per-SNP Fst values, set [window_size_in_bp] = 1 and [size_window_step_in_bp] = 1

To generate a circlize plot of the Fst values:
R package needed: <a href="https://cran.r-project.org/web/packages/circlize/index.html">circlize</a> (more details: <a href="https://jokergoo.github.io/circlize_book/book/">here</a>) 
<em>see ./3.2.7/Fst/Rscript_plot_circlize/script_circlize_Fst_africanrice.R</em>
</pre></code>

<pre><code>
<strong>Perform outlier detection (./3.2.7/pcadapt/Rscript_pcadapt/script_pcadapt.R) </strong>
R package needed: <a href="https://cran.r-project.org/web/packages/pcadapt/index.html">pcadapt</a> (more details: <a href="https://cran.r-project.org/web/packages/pcadapt/vignettes/pcadapt.html">here</a>) 
<em>To import the data, compute PCA, perform scans & generate Manhattan plots, see ./3.2.7/pcadapt/Rscript_pcadapt/script_pcadapt.R</em>
</pre></code>

### Pool-seq data - Application to the sessile oak data </br>
*3.3.1 - Pool-seq vs. individual data* <br/>
<pre><code>
Excel file needed: <a href="http://www1.montpellier.inra.fr/CBGP/software/PoolSeqUtils/">PIFs</a>  
<em>A simulation comparing the precision in allele frequency estimation for two sequencing strategies: pool-seq & individual sequencing (see <a href="https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.12360">here</a> for details).
Here a strategy based on a growing number of individuals sequenced at a pool coverage of 100X is compared to a strategy assuming 20 individuals sequenced separately at 20X. 
Results are shown in ./3.3.1/PIFs_simulation/Poolseq_vs_individual_100xpoolvs20xNbindividuals.sed. 
See also the R script ./3.3.1/PIFs_simulation/script_R_comp_power_poolseq_individual_sequencing.R. </em>
</pre></code>
*3.3.3 - From raw data to allele counts* <br/>
<pre><code>
<strong>1/ Downloading sequencing reads, trimming & reference genome indexing (see section 3.2.2 above)</strong>
Softwares needed: <a href="https://www.gnu.org/software/wget/">wget</a> (ftp-transfert) & <a href="https://github.com/timflutre/trimmomatic">Trimmomatic</a> (read trimming)
See <a href="https://www.ebi.ac.uk/ena/data/view/PRJEB21312">here</a> for a list of all sequencing data available (1 run accession = 1 lane, 4 lanes/pool). A correspondence table between SRA accession IDs and Population IDs as indicated in Leroy et al. 2019 is available (./3.3.3/IDs_correspondence_table.txt/)
The oak reference genome (PM1N) can be downloaded from <a href="http://www.oakgenome.fr/?page_id=587">oakgenome</a>.

<strong>2/ Mapping, sorting & removing duplicates (see ./3.3.3/4-Mapping/mapping_O16.sh)</strong>
Softwares needed: <a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">bowtie2</a>, <a href="http://samtools.sourceforge.net/">Samtools</a> (filtering & merging bams), <a href="https://broadinstitute.github.io/picard/">Picard</a> (sorting & removing duplicates)
To run bowtie2 with the sensitive end-to-end mode (see <a href="http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml">here</a> for more details): 
<em>bowtie2 -p [nb_CPUs_to_use] -k [search_alignments] -q --sensitive -x $ref -1 [file_containing_reads_end1] -2 [file_containing_reads_end2] -t --un-gz [output_file_reads_not_aligned] | samtools view -Shb /dev/stdin > [outfile.bam]
picard-tools/SortSam.jar INPUT= [outfile.bam] OUTPUT=[outfile.bam].pisorted SORT_ORDER=coordinate
picard-tools/MarkDuplicates.jar INPUT=[outfile.bam].pisorted OUTPUT=[outfile.bam].pisorted.dedup METRICS_FILE=[output_file_metrics] VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 REMOVE_DUPLICATES=true</em>
samtools merge [total_bam_pop1] [bam_lane1_pop1] [bam_lane2_pop1] [bam_lane3_pop1] [bam_lane4_pop1] 

<strong>3/ Pileup & synchronized pileup  (see ./3.3.3/5-Pileup_sync/script_samtools_mpileup.sh)</strong>
Softwares needed: <a href="http://samtools.sourceforge.net/">Samtools</a> & <a href="https://sourceforge.net/projects/popoolation2/">Popoolation2</a>
<em>samtools mpileup -f [fasta_genome] [total_bam_pop1] [total_bam_pop2] [total_bam_pop3] ... > [outfile_pileup]
java -Xmx4g -jar mpileup2sync.jar --input [outfile_pileup] --output [outfile_pileup].sync --fastq-type sanger --min-qual [minimum_base_quality] --threads [nb_CPUs_to_use]</em>

<strong>4/ Generate allele count for each position and SNP</strong>
Software needed: <a href="https://sourceforge.net/projects/popoolation2/">Popoolation2</a>
<em>popoolation/snp-frequency-diff.pl --input [outfile_pileup].sync --output-prefix [prefix_for_output_files] --min-count [minimum_number_of_alt_alleles] --min-coverage [min_coverage_per_pool] --max-coverage [max_coverage_per_pool]
./allele_count_withMAF.sh [previous_popoolation_prefix_for_output_files].rc [output-prefix] [tmpfile]</em>
This second script is just to give an example. You can also the R script <a href="https://cran.r-project.org/web/packages/poolfstat/index.html">poolfstat</a> to generate allele counts from the synchronized file, as described in the Rscript ./3.3.6/script_baypass/generate_baypass_inputfile.R (see section 3.3.6 below)
</pre></code>

*3.3.4 - Population splits & mixtures* <br/>
<pre><code>
<strong>Perform TreeMix simulations (./3.3.4/Script_TreeMix/treemix.sh)</strong>
Software needed: <a href="https://bitbucket.org/nygcresearch/treemix/wiki/Home">TreeMix</a>
<em>for i in {1..10}; do 
	./treemix-1.13/src/treemix -i [infile] -k 1000 -m $i -o [outfile].m.$i
done</em>

<strong>Compute explained variance, draw phylogenetic trees & show residuals (./3.3.4/Script_TreeMix/treemix.sh)</strong>
R scripts needed: <a href="https://bitbucket.org/nygcresearch/treemix/wiki/Home">TreeMix-associated R scripts</a> & <a href="https://cran.r-project.org/web/packages/ape/index.html">ape</a>
<em>See ./3.3.4/Rscript_TreeMix/script_R_treemix.R for details.</em>
</pre></code>

*3.3.5 - Fst estimates* <br/>
<pre><code>
<strong>Compute fixation index (Fst) from synchronized pileup using poolfstat(./3.3.5/)</strong>
R scripts needed: <a href="https://cran.r-project.org/web/packages/poolfstat/index.html">poolfstat</a>, <a href="https://cran.r-project.org/web/packages/reshape2/index.html">reshape2</a> & <a href="https://cran.r-project.org/web/packages/ggplot2/index.html">ggplot2</a>
<em>To import data from the popoolation2 synchronized mpileup format, compute pairwise population population Fst matrix, per-SNP pairwise or among-population Fst values or generate plots, see ./3.3.5/script_poolfstat_Fst_Hivert.R for details.</em>
</pre></code>
*3.3.6 - Genome Scan for Selection* <br/>

<pre><code>
<strong>Generate an infile (./3.3.6/script_baypass/script_baypass.sh)</strong>
R scripts needed: <a href="https://cran.r-project.org/web/packages/poolfstat/index.html">poolfstat</a>
From the popoolation2 synchronized mpileup format, the popsync2pooldata & pooldata2genobaypass functions are probably the easiest way to generate the baypass input file<em> (see ./3.3.6/script_baypass/generate_baypass_inputfile.R)</em>

<strong>Perform genome scans to detect candidate SNPs under selection (./3.3.6/script_baypass/script_baypass.sh)</strong>
Software needed: <a href="http://www1.montpellier.inra.fr/CBGP/software/baypass/">BayPass</a>
<em>i_baypass -npop [nb_pops] -gfile [INFILE] -poolsizefile [FILE_with_popsizes] [+ parameters related to Markov chains, read the BayPass's manual <a href="<a href="http://www1.montpellier.inra.fr/CBGP/software/baypass/files/BayPass_manual_2.1.pdf">here</a>">here</a>]</em>

<strong>Identify outlier loci (./3.3.6/Rscript_XtX/script_baypass_XtX_plots.R</strong>
R scripts needed: <a href="https://cran.r-project.org/web/packages/ggplot2/index.html">ggplot2</a>
<em>See ./3.3.6/Rscript_XtX/script_baypass_XtX_plots.R</em>
Generate pseudo-observed datasets (PODS) and Manhattan plots of XtX values highlighting outliers 
For more information on how to perform neutral simulations to calibrate the XtX, see also <a href="<a href="http://www1.montpellier.inra.fr/CBGP/software/baypass/files/BayPass_manual_2.1.pdf">here</a>">here</a>]</em>

</pre></code>
*3.3.7 - Genotype-Environment associations* <br/>
<pre><code>
<strong>Detect SNPs with clinal variation along environment or phenotypic gradients (./3.3.7/script_baypass_GEA-GPA.sh)</strong>
Software needed: <a href="http://www1.montpellier.inra.fr/CBGP/software/baypass/">BayPass</a>
<em>i_baypass -npop [nb_pops] -gfile [INFILE] -efile [FILE_with_evironmental] -poolsizefile [FILE_with_popsizes] -scalecov [+ parameters related to Markov chains, read the BayPass's manual <a href="<a href="http://www1.montpellier.inra.fr/CBGP/software/baypass/files/BayPass_manual_2.1.pdf">here</a>">here</a>]
It is important to use the "-scalecov" option to scale all covariables.</em>


<strong>Identify Genotype-environment (GEA) or Genotype-Phenotype associations (GPA/pGWAS) (./3.3.7/script_baypass_XtX_BFplots.R)</strong>
R scripts needed: <a href="https://cran.r-project.org/web/packages/ggplot2/index.html">ggplot2</a>
<em>See ./3.3.7/script_baypass_XtX_BFplots.R</em>
This script generates plots of XtX vs. Bayes Factors (BF) and Manhattan plots of BF values highlighting SNPs with strong support for associations with environmental or phenotypic variables
</pre></code>
