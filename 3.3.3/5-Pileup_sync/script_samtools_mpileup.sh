# TL - 180615
# generate a mpileup file using all merged files (see "mapping_O16.sh" for an example)
#$ -m abe
#$ -q unlimitq
cd /home/tleroy/work2/tleroy/
# samtools mpileup -f [reference file] [merged_bam_1] [merged_bam_2] [merged_bam_3] ... [merged_bam_18]   
samtools mpileup -f /home/tleroy/work2/tleroy/ref/haplome_v2.3.fa ./projet_BNP/A/BNP_AOSW.v2haplo.bam.pisorted.dedup.merged ./projet_BNP/B/BNP_BOSW.v2haplo.bam.pisorted.dedup.merged ./projet_BNP/C/BNP_COSW.v2haplo.bam.pisorted.dedup.merged ./projet_BNP/D/BNP_DOSW.v2haplo.bam.pisorted.dedup.merged ./projet_BNP/E/BNP_EOSW.v2haplo.bam.pisorted.dedup.merged ./projet_BNP/F/BNP_FOSW.v2haplo.bam.pisorted.dedup.merged ./projet_BNP/G/BNP_GOSW.v2haplo.bam.pisorted.dedup.merged ./projet_BNP/H/BNP_HOSW.v2haplo.bam.pisorted.dedup.merged ./projet_BNP/I/BNP_IOSW.v2haplo.bam.pisorted.dedup.merged ./projet_BNP/K/BNP_KOSW.v2haplo.bam.pisorted.dedup.merged ./projet_BHW/L/BHW_LOSW.v2.3.bam.pisorted.dedup.merged ./projet_BHW/I/BHW_IOSW.v2.3.bam.pisorted.dedup.merged ./projet_BHW/P/BHW_POSW.v2.3.bam.pisorted.dedup.merged ./projet_BHW/LH/BHW_LHOSW.v2.3.bam.pisorted.dedup.merged ./projet_BHW/GA/BHW_GAOSW.v2.3.bam.pisorted.dedup.merged ./projet_BHW/GE/BHW_GEOSW.v2.3.bam.pisorted.dedup.merged ./projet_BHW/A/BHW_AOSW.v2.3.bam.pisorted.dedup.merged ./projet_BHW/PR/BHW_PROSW.v2.3.bam.pisorted.dedup.merged > /home/tleroy/work2/tleroy/projet_BNP-BHW/BNP-BHW_18pops_v2.3.pileup


# Generate a synchronized mpileup (popoolation mpileup2syn)
java -Xmx4g -jar /usr/local/bioinfo/src/popoolation2/popoolation2_1201/mpileup2sync.jar --input /home/tleroy/work2/projet_BNP-BHW/BNP-BHW_18pops_v2.3.pileup --output /home/tleroy/work2/projet_BNP-BHW/BNP-BNP-BHW_18pops_v2.3.pileup.sync --fastq-type sanger --min-qual 20 --threads 24
