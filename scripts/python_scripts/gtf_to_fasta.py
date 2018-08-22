import sys
import os
import re
import copy
import argparse
import math
import gff_utils as gu
import bed_utils as bu
import fasta_utils as fu
from collections import Counter

########################
### ARGUMENT PARSING ###
########################

parser = argparse.ArgumentParser()
parser.add_argument('-A','--annotation',dest='reference_GFF',
                    help="Path to a reference annotation GTF/GFF.",required=True)
parser.add_argument('-G','--genome',dest='genome_fasta',
                    help="Path to a FASTA file of the genome.",required=True)
args = parser.parse_args()

#################
### FUNCTIONS ###
#################

def which(x,value=True):
    """Returns a list of locations in x that satisfy value"""
    return [a for a,b in enumerate(x) if b==value]


def notwhich(x,value=0):
    """Returns a list of locations in x that do not satisty value"""
    return [a for a,b in enumerate(x) if b!=value]


def flatten(list_of_lists):
    """Collapses a list/tuple of lists into a single list"""
    return [item for sublist in list_of_lists for item in sublist]

def std(list):
    """Calculates the standard deviation of a list of numbers"""
    if len(list) <= 1:
        return 0
    
    list = [float(i) for i in list]
    m = sum(list)/len(list)
    diffsq = [(i - m)**2 for i in list]
    return math.sqrt(sum(diffsq)/(len(list) - 1))


def is_inside(range_a,range_b):
    """Determines if one start/end double is contained
    within the range of the second"""
    return range_a[0] >= range_b[0] and range_a[1] <= range_b[1]


def convert_to_gene(dictionary, sep='.'):
    """ Collapse the values in a dictionary
    based on the key's prefix """
    keys = dictionary.keys()
    outdict = {}
    prefixes = list(set([k.split(sep)[0] for k in keys]))
    for i in range(len(prefixes)):
        combine = [k for k in keys if prefixes[i] in k]
        cval = 0
        for k in combine:
            cval += dictionary[k]
        outdict[prefixes[i]] = cval
    return outdict


#############################
### IMPORT REFERENCE DATA ###
#############################

strand_swap = {
    '+':'-',
    '-':'+',
    '.':'NA'
}

# Loads a FASTA file of the reference genome
genome = fu.import_genome(args.genome_fasta)
# print('# Genome FASTA: {}'.format(args.genome_fasta))

# Records the length of each chromosome in the reference
chromosomes = {}
for k,v in genome.items():
    chromosomes[k] = len(v)


if args.reference_GFF.split('.')[-1].lower() in ['gff','gtf','gff3']:
    ref_transcripts = gu.parse_annotation(args.reference_GFF)
elif args.reference_GFF.split('.')[-1].lower() in ['bed','bed12']:
    ref_transcripts = bu.read_bed12(args.ref_transcripts, 'BED12', 'reference', mode='dictionary')
        
ref_IDs = sorted(list(ref_transcripts.keys()))
# print('# {} reference transcripts: {}'.format(
    # len(ref_IDs),
    # args.reference_GFF
# ))

#################################
### EVALUATE THE GIVEN SAMPLE ###
#################################

for ID in ref_IDs:
    input_transcript = ref_transcripts[ID]
    chrom = input_transcript.chrom
    strand = input_transcript.strand
    # Gather information about exons
    exon_starts = input_transcript.get_exon_start()
    exon_ends = input_transcript.get_exon_end()
    if len(exon_starts) != len(exon_ends):
        print('# WARNING: inconsistent exon structure with {}'.format(ID))
        continue

    # Make a list of all nucleotide positions in the input transcript
    positions = flatten(
        [list(range(a - 1,b)) for a,b in zip(exon_starts,exon_ends)]
    )
    length = len(positions)
    positions = [i for i in positions if i < chromosomes[chrom] and i > 0]

    # Get transcript's nucleotide sequence from the FASTA file
    sequence = ''.join([genome[chrom][i] for i in positions])
    if strand == '-':
        sequence = fu.rc(sequence)
    
    print(
        '>{}\n{}'.format(
            ID,
            sequence
        )
    )
