# Introduction to population genomics
This github directory hosts all scripts used to perform the analyses shown in in the book chapter of Leroy & Rougemont (2019).

>*Leroy, T & Rougemont, Q. 2019. Population genetics analysis methods: from population structure inferences to genome scans for selection. In: Molecular Plant Taxonomy (Springer)*

Please send an email to both thibault.leroy@umontpellier.fr and quentinrougemont@orange.fr for questions regarding these scripts or for full-text requests.

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
<strong>Import sequencing data (./1-Import_RawData0)</strong>
Softwares needed: wget (ftp-transfert)
<em> <a href="https://www.ebi.ac.uk/ena/data/view/PRJEB21312"> Have a look here </a>
e.g.wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR200/000/ERR2008850/ERR2008850_1.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR200/000/ERR2008850/ERR2008850_2.fastq.gz</em>
</code></pre>



<pre><code>
<strong>Read trimming (./2-Trimming)</strong>
Softwares needed: Trimmomatic
<em>java -Xmx4g -jar ./trimmomatic-0.33.jar PE -threads 1 -phred33 "[file]_1.fastq.gz [file]_2.fastq.gz [file]_1_cleaned.fastq.gz [file]_1_cleaned_unpaired.fastq.gz [file]_2_cleaned.fastq.gz [file]_2_cleaned_unpaired.fastq.gz ILLUMINACLIP:./adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50</em>
</pre></code>

- 
