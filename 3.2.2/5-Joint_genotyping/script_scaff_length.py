#!/usr/bin/python
# TL - 160218
from Bio import SeqIO
import sys
fastafile=str(sys.argv[1])
for seq_record in SeqIO.parse(fastafile, "fasta"):
 printed_line = '%s\t%i' % (seq_record.id, len(seq_record))
 print(printed_line)
