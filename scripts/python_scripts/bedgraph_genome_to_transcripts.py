import re
import sys
import argparse
import fasta_utils as fu
import numpy as np

########################
### ARGUMENT PARSING ###
########################

desc = (
    "Takes a genome-mapped BEDGRAPH file and an annotation file (GFF),"
    "outputs a new BEDGRAPH file where 'chromosomes' are each transcript."
)

parser = argparse.ArgumentParser()
parser.add_argument("bedgraph_plus",
                    help="BEDGRAPH of plus-stranded nucleotide-level coverage.")
parser.add_argument("bedgraph_minus",
                    help="BEDGRAPH of minus-stranded nucleotide-level coverage.")
parser.add_argument(dest='reference_GFF',
                    help="Path to a reference annotation GTF/GFF.")
parser.add_argument(dest='genome_fasta',
                    help="Path to a FASTA file of the genome.")
parser.add_argument('--subset', default='', type=str,
                    help="Path to file of transcript IDs to use from reference.")
parser.add_argument('-O','--output', default='transcript_coverage.bedgraph', type=str,
                    help="Filename to output the transcript-level bedgraph.")
parser.add_argument('--write_fasta', default='', type=str,
                    help="(optional) Filename to output the transcript FASTA file.")
parser.add_argument('--delimiter', default='\.', type=str,
                    help="Character(s) delimiting gene name from isoform number.")

args = parser.parse_args()

####################
# DEFINE FUNCTIONS #
####################

def which(x,value=True):
    """Returns a list of locations in x that satisfy value"""
    return [a for a,b in enumerate(x) if b==value]


def flatten(list_of_lists):
    """Collapses a list/tuple of lists into a single list"""
    return [item for sublist in list_of_lists for item in sublist]

def vector_to_bedgraph_lines(value_vector, chrom, chromstart=0, digits=2, multi_key=[]):
    """Converts a vector of values to bedgraph format, yielding lines."""
    start = 0
    prevpos = 0
    prevcount = None
    position = None
    vpos = None
    # Iterate through every position with values in the dict
    
    for vpos in [k for k,v in enumerate(value_vector) if v != 0]:
        if multi_key:
            all_values = [str(round(value_vector[vpos].get(k,0),digits)) for k in multi_key]
            count = '\t'.join(all_values)
        else:
            count = round(value_vector[vpos],digits)
        
        if count != prevcount or int(vpos) > 1 + prevpos:
            # The newly encountered value is not a continuation
            # of the previous value. Write the old run and start another.
            if prevcount and prevcount != 0:
                line_to_write = '\t'.join(
                    [
                        str(i) for i in [
                            chrom,
                            start + chromstart,
                            prevpos + 1 + chromstart,
                            prevcount
                        ]
                    ]
                )
                
                yield line_to_write
            
            start = vpos
        
        prevcount = count
        prevpos = int(vpos)

    if vpos and prevcount and prevcount != 0:
        line_to_write = '\t'.join(
            [
                str(i) for i in [
                    chrom,
                    start + chromstart,
                    prevpos + 1 + chromstart,
                    prevcount
                ]
            ]
        )
        
        yield line_to_write

####################
# LOAD ENVIRONMENT #
####################

# TODO: retool gff_utils to import gff3 format files
ref_transcripts = {}
refGFF = open(args.reference_GFF)
delim = args.delimiter
for line in refGFF:
    if line[0] == '#':continue
    chrom,source,gtype,start,end,score,strand,phase,other = line.rstrip().split('\t')
    if gtype not in ['exon','pseudogenic_exon','miRNA_primary_transcript']:
        continue
    
    parent = re.search('^.*;?Parent=([^;]+);?.*$',other).groups()[0]
    ID = re.search('(gene|transcript)_id=(.+?)[\.:;].+$',other)
    if ID:
        ID = ID.groups()[0]
    else:
        ID = parent.split(delim)[0]
    
    iso_length = len(re.findall(',',parent))+1
    isos = re.search(','.join([ID+'('+delim+'[0-9]+)?']*iso_length),parent).groups()
    iso_ID = parent.split(',')
    for i in iso_ID:
        ref_transcripts[i] = ref_transcripts.get(i,{})
        if gtype != 'miRNA_primary_transcript':
            exon_number = re.search(
                'ID=' + ID + ':(pseudogenic_)?exon:([0-9]+).*;Parent.+$',
                other
            )
            if exon_number:
                exon_number = exon_number.groups()[1]
            else:
                exon_number = str(max([int(j) for j in ref_transcripts[i].get('exon_nums',['0'])])+1)
        else:
            exon_number = '1'
    
    for i in isos:
        if i:
            iso_ID = ID + i
        else:
            iso_ID = ID
        
        ref_transcripts[iso_ID]['chrom'] = chrom
        ref_transcripts[iso_ID]['strand'] = strand
        ref_transcripts[iso_ID]['start'] = ref_transcripts[iso_ID].get('start',[])+[int(start)]
        ref_transcripts[iso_ID]['end'] = ref_transcripts[iso_ID].get('end',[])+[int(end)]
        ref_transcripts[iso_ID]['exon_nums'] = ref_transcripts[iso_ID].get('exon_nums',[])+[exon_number]  
        
ref_IDs = sorted(list(ref_transcripts.keys()))
print('# {} reference transcripts: {}'.format(
    len(ref_IDs),
    args.reference_GFF
))

# 'picked_IDs' is an array of IDs to use from the reference_GFF
if args.subset:
    picked_IDs = [i.rstrip().split('\t')[0] for i in open(args.subset).readlines()]
else:
    picked_IDs = ref_IDs

# 'genome' is a dict of strings for each chromosome in 'genome_fasta'
genome = fu.import_genome(args.genome_fasta)

# 'chromosomes' contains the lengths of all chromosomes the that BEDGRAPH contains values for.
chromosomes = {}
for chrom in genome.keys():
    length = len(genome[chrom])
    chromosomes[chrom] = int(length)

# 'coverage' is a dictionary of float vectors for each nucleotide in the genome.
# Contains the value of the BEDGRAPH file at each position.
coverage = {}
coverage['+'] = {}
coverage['-'] = {}

for chrom,chromlen in chromosomes.items():
    coverage['+'][chrom] = np.zeros(chromlen, dtype='float32')
    coverage['-'][chrom] = np.zeros(chromlen, dtype='float32')

# Populate the + bedgraph dict
coverage_file = open(args.bedgraph_plus)
for line in coverage_file:
    chrom,start,end,count = line.rstrip().split()
    count = float(count)
    coverage['+'][chrom][int(start):int(end)] += count

# Populate the - bedgraph dict
coverage_file = open(args.bedgraph_minus)
for line in coverage_file:
    chrom,start,end,count = line.rstrip().split()
    count = float(count)
    coverage['-'][chrom][int(start):int(end)] += count

# Open a FASTA file for the picked IDs
if args.write_fasta:
    fasta_outfile = open(args.write_fasta,'w')
else:
    del genome

# Generate an output dict of transcript-level coverages
output_file = open(args.output, 'w')

print('# Converting genome BEDGRAPH to transcript coordinates.')
for ID in picked_IDs:
    if ID not in ref_transcripts:
        print('# WARNING: {} not found'.format(ID))
        continue
    
    chrom = ref_transcripts[ID]['chrom']
    strand = ref_transcripts[ID]['strand']
    exon_starts = sorted(ref_transcripts[ID]['start'])
    exon_ends = sorted(ref_transcripts[ID]['end'])

    transcript_coverage = np.empty(0, dtype='float32')
    if args.write_fasta:
        transcript_fasta = ''

    if strand not in ['+','-']:
        print('# WARNING: {} is unstranded'.format(ID))
        continue

    for a,b in zip(exon_starts, exon_ends):
        exon = coverage[strand][chrom][a-1:b]
        transcript_coverage = np.append(transcript_coverage, exon)
        if args.write_fasta:
            transcript_fasta += genome[chrom][a-1:b]
    
    
    if strand == '-': # Flip the strand for minus-stranded data
        transcript_coverage = transcript_coverage[::-1]
        if args.write_fasta:
            transcript_fasta = fu.rc(transcript_fasta)

    for line in vector_to_bedgraph_lines(transcript_coverage, ID):
        output_file.write(line+'\n')
    
    if args.write_fasta:
        fasta_outfile.write('>{}\n'.format(ID))
        fasta_outfile.write('{}\n'.format(transcript_fasta))

output_file.close()
if args.write_fasta:
    fasta_outfile.close()
