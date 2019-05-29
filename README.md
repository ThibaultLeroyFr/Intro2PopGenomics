# Introduction to population genomics
This github directory hosts all scripts used to perform the analyses shown in in the book chapter of Leroy & Rougemont (2019).

>*Leroy, T & Rougemont, Q. 2019. Population genetics analysis methods: from population structure inferences to genome scans for selection. In: Molecular Plant Taxonomy (Springer)*

Please send an email to <a href="mailto:thibault.leroy@umontpellier.fr;quentinrougemont@orange.fr?subject=[Intro2popGenomics-Github]">us </a> and for questions regarding these scripts or for full-text requests.

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
Our scripts are not standalone executables. Quite the contrary, these scripts (deliberately) require some simple edits to adjust to your data. Editing script is probably the best way to learn how a script works, to detect and correct the errors and, more broadly, to start learning how to code. So please take some time to read the scripts and, ideally, to briefly look at the software user manuals. Change the paths to files and programs to adjust the scripts to your computer architecture. Do not hesitate to send emails in case of major computational issues. 


## Details
### Individual data:</br>
*3.2.2 : From raw data to VCF*</br>

<pre><code>
<strong>1/Import sequencing data (./3.2.2/1-Import_RawData/)</strong>
Softwares needed: <a href="https://www.gnu.org/software/wget/">wget</a> (ftp-transfert)
<em><a href="https://www.ebi.ac.uk/ena/data/view/PRJEB21312">Have a look here to get a list of the ftp files</a>
e.g.wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR200/000/ERR2008850/ERR2008850_1.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR200/000/ERR2008850/ERR2008850_2.fastq.gz</em>
</code></pre>

<pre><code>
<strong>2/Read trimming (./3.2.2/2-Trimming/)</strong>
Softwares needed: <a href="https://github.com/timflutre/trimmomatic">Trimmomatic</a>
<em>java -Xmx4g -jar ./trimmomatic-0.33.jar PE -threads 1 -phred33 "[file]_1.fastq.gz [file]_2.fastq.gz [file]_1_cleaned.fastq.gz [file]_1_cleaned_unpaired.fastq.gz [file]_2_cleaned.fastq.gz [file]_2_cleaned_unpaired.fastq.gz ILLUMINACLIP:./adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50</em>
</pre></code>

<pre><code>
<strong>3/Downloading & indexing a reference genome (./3.2.2/3-Download_Index_references/)</strong>
Softwares needed: <a href="https://www.gnu.org/software/wget/">wget</a> (ftp-transfert), <a href="https://sourceforge.net/projects/bio-bwa/files/">BWA</a>, <a href="http://samtools.sourceforge.net/">Samtools</a> & <a href="https://broadinstitute.github.io/picard/">Picard</a> (indexing) 
<em>Asian rice genome: wget ftp://ftp.ensemblgenomes.org/pub/plants/release-42/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna_sm.toplevel.fa.gz
BWA: bwa index -a bwtsw [fasta]
Samtools: samtools faidx [fasta]
Picard: java -Xmx4g  -jar picard.jar  CreateSequenceDictionary REFERENCE=[fasta]  OUTPUT=[fasta].dict </em>
</pre></code>

<pre><code>
<strong>4/Mapping & Individual Calling (./3.2.2/4-PipelineMappingCalling/)</strong>
Softwares needed: <a href="https://sourceforge.net/projects/bio-bwa/files/">BWA</a> (mapping), <a href="http://samtools.sourceforge.net/">Samtools</a> (filtering), <a href="https://broadinstitute.github.io/picard/">Picard</a> (removing duplicates) & <a href="https://software.broadinstitute.org/gatk/download/">GATK</a> (creating gVCF files)
<em>bash 1_mapping.sh [RawData_Directory] [File_Trim_Paired_1] [File_Trim_Paired_2] [File_Trim_Unpaired_1] [File_Trim_Unpaired_2] [Reference_Genome] [Number_of_CPU_to_use]
bash 2_snpindel_callingGVCF.sh [ID] [Reference_Genome] [output_directory] [Number_of_CPU_to_use] </em>
</pre></code>

<pre><code>
<strong>5/Joint Genotyping (./3.2.2/5-Joint_genotyping/)</strong>
<em> bash 3_intervals_jointgenotyping.sh [reference_genome] [nb_cpus]  </em>
</pre></code>

<pre><code>
<strong>6/Filtering variants (./3.2.2/6-Variant_Filtration/)</strong>
<em> ./VariantFiltrationVCF.py -q 2.0 -s 60.0 -m 40.0 -n -2.0 -r -2.0 -w 45000 -f [VCF] > [filtered_VCF] </em>
</pre></code>

*3.2.3 - Population structure*</br>
<pre><code>
<strong>Perform Principal Component Analysis (./3.2.3/PCA/)</strong>
R packages needed: <a href="https://bioconductor.org/packages/release/bioc/html/SNPRelate.html/">SNPRelate</a> (& <a href="https://bioconductor.org/packages/release/bioc/html/gdsfmt.html/">gdsfmt</a>)
<em> Details shown in script_PCA_from_vcf.R (Import & convert the vcf file, compute PCA & generate plots) </em>
</pre></code>

<pre><code>
<strong>Infer ancestry proportions (./3.2.3/fastStructure/)</strong>
Software needed: <a href="https://www.cog-genomics.org/plink2/">Plink</a> & <a href="https://rajanil.github.io/fastStructure/">fastStructure</a>
Generate an input file using Plink:
<em> ./plink --allow-extra-chr -allow-extra-chr --make-bed --noweb --out [VCF].bed --vcf [VCF] </em>
FastStructure inferences:
<em>./structure.py -K [nb_of_clusters] --input [VCF].bed --output [VCF].bed.K[nb_of_clusters] --full --cv [cross-validation steps] --format bed
./chooseK.py --input=[VCF].bed.K
</pre></code>

*3.2.4 - Nucleotide diversity*<br/>
<pre><code>
<strong> Generate Fasta files from VCF (./3.2.4/Scripts_generate_fasta_sequences_from_vcf/)</strong>
Scripts were developed by Leroy et al. 2019 bioRxiv 505610, ver. 4 peer-reviewed and recommended by PCI Evolutionary Biology (paper: <a href="https://www.biorxiv.org/content/biorxiv/early/2019/05/24/505610.full.pdf">here</a>, scripts: <a href="https://figshare.com/s/122efbec2e3632188674">here</a>)
<em>python ../VCF2Fasta_fast.py -q [threshold_base_quality] -m [minimum_coverage_position] -M [maximum_coverage_position] -f [VCF] > [outputdir]/Outputs_VCF2Fasta.txt</em>

<strong> Compute theta, pi & Tajima's D (./3.2.4/Scripts_generate_fasta_sequences_from_vcf/)</strong>
Software needed: <a href="https://figshare.com/s/122efbec2e3632188674#/articles/7484705">seq_stat</a>
<em>bash script_compute_pi_slidwin.sh [referencegenome].scafflength [VCF] [ouputdir] [size_of_sliding_window] [output_prefix]</em>
</pre></code>
