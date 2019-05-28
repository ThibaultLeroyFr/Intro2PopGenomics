#!/usr/bin/env python3
"""Extract genomic regions around SNPs and create 2 variants

Usage:
    <program> input_genome wanted_file flanking_size output_fasta

Where:
    input_genome is a fasta or fasta.gz genome file
    wanted_file has this format:
    - Has 5 tab-separated columns
        - Scaffold name (as in genome file)
        - Position on scaffold
        - SNPID_POS
        - SNP Allele 1 (major)
        - SNP Allele 2 (minor)
    - May contain header files starting with '#'

    #CHROM	POS	ID	REF	ALT
    gi|1134601138|ref|NC_025992.3|	17342	664_11	A	C
    gi|1134601138|ref|NC_025992.3|	17346	664_15	T	A
    gi|1134601138|ref|NC_025992.3|	17361	664_30	A	C
    gi|1134601138|ref|NC_025992.3|	17378	664_47	C	T
    gi|1134601138|ref|NC_025992.3|	17387	664_56	C	A
    gi|1134601138|ref|NC_025992.3|	17388	664_57	A	T

    flanking_size is the number (int) of nucleotides to keep on each side
    output_fasta is the name of the extracted SNP variant fasta file
"""

# Modules
from collections import defaultdict
import gzip
import sys

# Defining classes
class Fasta(object):
    """Fasta object with name and sequence
    """

    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

    def write_to_file(self, handle):
        handle.write(">" + self.name + "\n")
        handle.write(self.sequence + "\n")

    def __repr__(self):
        return self.name + " " + self.sequence[0:31]

# Functions
def fasta_iterator(input_file):
    """Takes a fasta file input_file and returns a fasta iterator
    """
    with myopen(input_file) as f:
        sequence = ""
        name = ""
        begun = False
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if begun:
                    yield Fasta(name, sequence)
                name = line.replace(">", "")
                sequence = ""
                begun = True
            else:
                sequence += line

        if name != "":
            yield Fasta(name, sequence)

def myopen(infile, mode="rt"):
    """Replacement for `open` function to accept gzip files

    Use gzip compression algorithm on files ending with `.gz`

    `myopen` can be used like the `open` function because it has the same
    behaviour. Namely, it returns a file handle (ie: an opened connection
    to a file).
    """

    # If filename ends with .gz, open in gzip mode
    if infile.endswith(".gz"):
        return gzip.open(infile, mode=mode)

    # Else open normally
    else:
        return open(infile, mode=mode)

# Parse user input
try:
    input_genome = sys.argv[1]
    wanted_file = sys.argv[2]
    flanking_size = int(sys.argv[3])
    output_fasta = sys.argv[4]
except:
    print(__doc__)
    sys.exit(1)

# Create wanted set
wanted_regions = defaultdict(list)

with open(wanted_file) as wfile:
    for line in wfile:
        # Skip comment lines
        if line.startswith("#"):
            continue

        scaffold, position, snp_id, allele1, allele2 = line.strip().split("\t")
        wanted_regions[scaffold].append((position, snp_id, allele1, allele2))

# Extract regions
sequences = fasta_iterator(input_genome)

with open(output_fasta, "w") as outfile:
    count_good_position = 0
    count_bad = 0

    for scaffold in sequences:
        scaffold.name = scaffold.name.split(" ")[0]
        scaffold_id = scaffold.name.split()[0]

        if scaffold_id in wanted_regions:
            seq = scaffold.sequence

            for snp in wanted_regions[scaffold_id]:

                # Compute SNP position
                position, snp_id, allele1, allele2 = snp

                pos = int(position) - 1

                assert pos <= len(seq), "Errof: SNP outside the scaffold"

                left = pos - flanking_size

                if left < 0:
                    left = 0

                right = pos + flanking_size + 1

                # Create variant sequences
                variant_seq1 = seq[left: pos].upper() + allele1 + seq[pos + 1: right]
                variant_seq2 = seq[left: pos].upper() + allele2 + seq[pos + 1: right]

                variant1 = Fasta("_".join([scaffold_id, position, snp_id, allele1]),
                        variant_seq1)

                variant2 = Fasta("_".join([scaffold_id, position, snp_id, allele2]),
                        variant_seq2)

                # Write to output file
                variant1.write_to_file(outfile)
                variant2.write_to_file(outfile)
