# Introduction to population genomics
This github directory hosts all scripts used to perform the analyses shown in in the book chapter of Leroy & Rougemont (2019).

>*Leroy, T & Rougemont, Q. 2019. Population genetics analysis methods: from population structure inferences to genome scans for selection. In: Molecular Plant Taxonomy (Springer)*

Please send an email to <a href="mailto:thibault.leroy@umontpellier.fr;quentinrougemont@orange.fr?subject=[Intro2popGenomics-Github]">us </a> for questions regarding these scripts or for full-text requests.

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


## Details
### Individual data:</br>
*3.2.2 : From raw data to VCF*</br>

<pre><code>
<strong>1/Import sequencing data (./3.2.2/1-Import_RawData/)</strong>
Softwares needed: <a href="https://www.gnu.org/software/wget/">wget</a> (ftp-transfert)
<em><a href="https://www.ebi.ac.uk/ena/data/view/PRJEB21312">Have a look here to get a list containing all ftp addresses</a>
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
(the file containing the summary statistics was made available at ./3.2.4/Results_pi/Genomic_pi_ROD_010419.withoutNA.txt)
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
<em>./script_parser_vcf.py [VCF_focal_species_ONLY_PASS_variants] [VCF outgroup1] [VCF outgroup2]> [Merged_VCF_file]<\em>

Then use awk '{print $X"    "Y...}' (where X and Y correspond to the columns in the [Merged_VCF_file]) to parse the data to obtain the following file format:
<em>chr    pos focal_All1   focal_All2  outgroup1_all1  outgroup1_all2 outgroup2_all1   outgroup2_all2
1	1248	G	A	G	.	G	.<\em>

Detect the ancestral allele & compute allele frequencies
<em>./script_ancestral_derived_counts.py [this_infile] > [file_with_derived_counts]
./script_compute_derivedallelefreq.py [file_with_derived_counts] > [files_with_derivedallfreq]<\em>

</pre></code>


3.2.7 - Fst & genome scans
<pre><code>

</pre></code>

### Pool-seq data:</br>
*3.3.1 - Pool-seq vs. individual data* <br/>
<pre><code>

</pre></code>
*3.3.3 - From raw data to allele counts* <br/>
<pre><code>

</pre></code>
*3.3.4 - Population splits & mixtures* <br/>
<pre><code>

</pre></code>
*3.3.5 - Fst estimates* <br/>
<pre><code>

</pre></code>
*3.3.6 - Genome Scan for Selection* <br/>
<pre><code>

</pre></code>
*3.3.7 - Genotype-Environment associations* <br/>
<pre><code>

</pre></code>
