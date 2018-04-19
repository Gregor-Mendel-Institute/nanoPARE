import sys
import os
import re
import math
import argparse
import fasta_utils as fu
from collections import Counter

desc = (
    "Takes a BEDGRAPH, as well as BED file(s) with upstream/downstream regions"
    "to set to zero; outputs masked BEDGRAPH."
    "Any values immediately upstream of -U or downstream of -D will be set to zero."
    "This can be used to mask sequence-specific artifacts, like TSO strand invasion or oligo-dT mispriming."
    "Any values inside features of the -I bed file will be set to zero."
)
parser = argparse.ArgumentParser(description=desc)

# add arguments to the ArgumentParser
parser.add_argument(
    '-L','--lengths', dest='lengths', type=str, 
    help='filepath to chromosome lengths table',
    required=True
)
parser.add_argument(
    '-I', '--input', dest='input', type=str, 
    help='input bedgraph coverage file(s)',
    default=[], nargs='+'
)
parser.add_argument(
    '-N', '--names', dest='names', type=str, 
    help='names corresponding to input bedgraphs',
    default=[], nargs='+'
)
parser.add_argument(
    '-F', '--features', dest='features', type=str, 
    help='input feature bed file',
    required=True
)
parser.add_argument(
    '-O', '--output', dest='output', type=str, 
    help='output table filepath',
    default='feature_counts.tsv'
)
parser.add_argument(
    '--bed_out', dest='BED_OUT', action='store_true', 
    help='output a BED format file',
    default=False
)
parser.add_argument(
    '-G', '--genome', dest='genome', type=str, 
    help='input genome fasta file',
    default=None
)
parser.add_argument(
    '--g_content', type=str, default=None,
    help='filepath to output G content table'
)

args = parser.parse_args()

def which(x,value=True):
    return [a for a,b in enumerate(x) if b==value]


#'chromosomes' contains the lengths of all chromosomes the that BEDGRAPH contains values for.
# Expects a two-column tab-separated file with:
#    chromosome  length
# Provided with the 'lengths' argument.
if args.genome:
    genome = fu.import_genome(args.genome)

chromosomes={}
lengths_file=open(args.lengths)
for line in lengths_file:
    chrom,length=line.rstrip().split('\t')
    chromosomes[chrom]=int(length)

#'coverage' is a nested dictionary of float vectors for each nucleotide in the genome.
# Contains a dictionary for each BEDGRAPH file with values at each position.
coverage = {}

ingraphs = args.input
for graph in ingraphs:
    print('Importing {}...'.format(graph))
    coverage[graph] = {}
    for chrom,chromlen in chromosomes.items():
        coverage[graph][chrom] = {}
    coverage_file = open(graph)
    for line in coverage_file:
        if line[0] == '#':
            continue
        
        chrom,start,end,count = line.rstrip().split()
        count = float(count)
        for i in range(int(start),int(end)):
            coverage[graph][chrom][i] = count

bedfile=open(args.features)
bed_features = {}
print('Loading bed features...')
strandcol = 5
for line in bedfile:
    if line[0]=='#':
        continue
    
    l = line.rstrip().split('\t')
    chrom = l[0]
    start = int(l[1])
    end = int(l[2])+1
    name = l[3]
    strand = l[strandcol]
    if name not in bed_features:
        bed_features[name] = dict(
            [
                ('chrom',chrom),
                ('strand',strand),
                ('positions',set(range(start,end)))
            ]
        )
    else:
        bed_features[name]['positions'].update(
            set(range(start,end))
        )

print('Calculating coverages...')
outfile = open(args.output, 'w')

outlengths = open('feature_lengths.tsv','w')
if args.g_content:
    g_file = open(args.g_content, 'w')

if not args.BED_OUT:
    outfile.write('feature\tchrom\tstrand\t'+'\t'.join(args.names)+'\n')

for feature in sorted(bed_features.keys()):
    outlengths.write(
        '{}\t{}\n'.format(
            feature,
            len(bed_features[feature]['positions'])
        )
    )
    
    values = [
        sum(
            [
                coverage[graph][bed_features[feature]['chrom']].get(i,0)
                for i in bed_features[feature]['positions']
            ]
        )
        for graph in ingraphs
    ]
    
    if args.g_content:
        nucleotides = [
            genome[bed_features[feature]['chrom']][i] 
            for i in bed_features[feature]['positions']
        ]
        if bed_features[feature]['strand'] == '-':
            nucleotides = [i for i in fu.rc(''.join(nucleotides))]
        
        G = round(float(sum([i=='G' for i in nucleotides])) / len(nucleotides),3)
        g_file.write('{}\t{}\n'.format(
            feature,
            G
        ))
    if args.BED_OUT:
        outfile.write(
            '\t'.join(
                [
                    str(i)
                    for i in [
                        bed_features[feature]['chrom'],
                        min(bed_features[feature]['positions']),
                        max(bed_features[feature]['positions']),
                        feature,
                        '.',
                        bed_features[feature]['strand'],
                    ] + values
                ]
            ) + '\n'
        )    
    else:
        outfile.write(
            '\t'.join(
                [
                    str(i)
                    for i in [
                        feature,
                        bed_features[feature]['chrom'],
                        bed_features[feature]['strand'],
                    ] + values
                ]
            ) + '\n'
        )

outfile.close()
outlengths.close()
if args.g_content:
    g_file.close()