import sys
import os
import re
import math
import argparse
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
    '-P', '--plus', dest='PLUS_BG', type=str, 
    help='input plus strand bedgraph filepath',
    default=None
)
parser.add_argument(
    '-M', '--minus', dest='MINUS_BG', type=str, 
    help='input minus strand bedgraph filepath',
    default=None
)
parser.add_argument(
    '-PO', '--plus_out', dest='PLUS_OUT', type=str, 
    help='output plus strand bedgraph filepath',
    default=None
)
parser.add_argument(
    '-MO', '--minus_out', dest='MINUS_OUT', type=str, 
    help='output minus strand bedgraph filepath',
    default=None
)
parser.add_argument(
    '-U', '--upstream', dest='BED_UP', type=str, 
    help='Bed file of features to mask upstream of',
    default=None
)
parser.add_argument(
    '-D', '--downstream', dest='BED_DOWN', type=str, 
    help='Bed file of features to mask downstream of',
    default=None
)
parser.add_argument(
    '-I', '--inside', dest='BED_INSIDE', type=str, 
    help='Bed file of features to mask inside of',
    default=None
)

args = parser.parse_args()

def which(x,value=True):
    return [a for a,b in enumerate(x) if b==value]

'''
'chromosomes' contains the lengths of all chromosomes the that BEDGRAPH contains values for.
Expects a two-column tab-separated file with:
    chromosome  length
Provided with the -L argument.
''' 
chromosomes={}
lengths_file=open(args.lengths)
for line in lengths_file:
    chrom,length=line.rstrip().split('\t')
    chromosomes[chrom]=int(length)

mask_positions={}
for chrom in list(chromosomes.keys()):
    mask_positions[chrom]={}
    mask_positions[chrom]['-']=set()
    mask_positions[chrom]['+']=set()

if args.BED_UP:
    print('Upstream mask: '+args.BED_UP)
    mask_file=open(args.BED_UP)    
    for line in mask_file:
        if line[0]=='#':continue
        l=line.rstrip().split()
        chrom=l[0]
        start_pos=int(l[1])
        end_pos=int(l[2])
        strand = l[5]
        if strand=='+':
            mask_positions[chrom][strand].add(end_pos)
        elif strand=='-':
            mask_positions[chrom][strand].add(start_pos-1)
    mask_file.close()

if args.BED_DOWN:
    print('Downstream mask: '+args.BED_DOWN)
    mask_file=open(args.BED_DOWN)    
    for line in mask_file:
        if line[0]=='#':continue
        l=line.rstrip().split()
        chrom=l[0]
        start_pos=int(l[1])
        end_pos=int(l[2])
        strand = l[5]
        if strand=='+':
            mask_positions[chrom][strand].add(start_pos-1)
        elif strand=='-':
            mask_positions[chrom][strand].add(end_pos)
    mask_file.close()
    
if args.BED_INSIDE:
    print('Internal mask: '+args.BED_INSIDE)
    mask_file=open(args.BED_INSIDE)    
    for line in mask_file:
        if line[0]=='#':continue
        l=line.rstrip().split()
        chrom=l[0]
        start_pos=int(l[1])
        end_pos=int(l[2])
        strand = l[5]
        mask_positions[chrom][strand].update(set(range(start_pos,end_pos)))
    mask_file.close()

if args.PLUS_BG:
    bedgraph_file=open(args.PLUS_BG)
    bedgraph_outfile=open(args.PLUS_OUT,'w')
    
    for line in bedgraph_file:
        all_ranges = None
        positions = None
        l = line.split('\t')
        chrom = l[0]
        start = int(l[1])
        end = int(l[2])
        if chrom not in chromosomes:
            continue
        
        positions = [
            i
            for i in range(start,end)
            if not i in mask_positions[chrom]['+']
        ]
        
        if not positions:
            continue
        
        if len(positions) == 1:
            all_ranges = [(start,end)]
        else:
            breaks = which(
                [
                    a - b != 1
                    for a,b in zip(
                        positions[1:],
                        positions[:-1]
                    )
                ]
            )
            if breaks:
                partitions = [0] + breaks + [len(positions)-1]
                all_ranges = [
                    (positions[a],positions[b]+1)
                    for a,b in zip(
                        partitions[:-1],
                        partitions[1:]
                    )
                ]
            else:
                all_ranges = [(start,end)]
        if all_ranges:
            for s,e in all_ranges:
                if len(l) > 3:
                    newline = '{}\t{}\t{}\t{}'.format(
                        chrom,
                        start,
                        end,
                        '\t'.join(l[3:])
                    )
                else:
                    newline = '{}\t{}\t{}\n'.format(
                        chrom,
                        start,
                        end
                    )
                bedgraph_outfile.write(newline)
            
    bedgraph_outfile.close()

if args.MINUS_BG:
    bedgraph_file=open(args.MINUS_BG)
    bedgraph_outfile=open(args.MINUS_OUT,'w')
    
    for line in bedgraph_file:
        all_ranges = None
        positions = None
        l = line.split('\t')
        chrom = l[0]
        start = int(l[1])
        end = int(l[2])
        if chrom not in chromosomes:
            continue
        
        positions = [
            i
            for i in range(start,end)
            if not i in mask_positions[chrom]['-']
        ]
        
        if not positions:
            continue
        
        if len(positions) == 1:
            all_ranges = [(start,end)]
        else:
            breaks = which(
                [
                    a - b != 1
                    for a,b in zip(
                        positions[1:],
                        positions[:-1]
                    )
                ]
            )
            if breaks:
                partitions = [0] + breaks + [len(positions)-1]
                all_ranges = [
                    (positions[a],positions[b]+1)
                    for a,b in zip(
                        partitions[:-1],
                        partitions[1:]
                    )
                ]
            else:
                all_ranges = [(start,end)]
        
        if all_ranges:
            for s,e in all_ranges:
                if len(l) > 3:
                    newline = '{}\t{}\t{}\t{}'.format(
                        chrom,
                        start,
                        end,
                        '\t'.join(l[3:])
                    )
                else:
                    newline = '{}\t{}\t{}\n'.format(
                        chrom,
                        start,
                        end
                    )
                bedgraph_outfile.write(newline)
    
    bedgraph_outfile.close()
