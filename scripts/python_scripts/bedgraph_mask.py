import sys,os,re,math
from collections import Counter

usage="\
Takes a BEDGRAPH, as well as BED file(s) with upstream/downstream regions to set to zero; outputs masked BEDGRAPH.\
\n\
Commandline arguments:\n\
    -B=[input bedgraph]           (required filepath)\n\
    -O=[output bedgraph]          (required filepath)\n\
    -L=[lengths table]            (required filepath; output from fasta_lengths.py)\n\
    -S=[strand]                   (default: both; options: both|. , plus|+ , minus|-)\n\
    -U=[upstream BED mask]        (filepath)\n\
    -D=[downstream BED mask]      (filepath)\n\
\n\
Any values immediately upstream of -U or downstream of -D will be set to zero.\n\
This can be used to mask sequence-specific artifacts, like TSO strand invasion or oligo-dT mispriming.\
"

BED_UP   = 'none'
BED_DOWN = 'none'
STRAND   = 'both'

if len(sys.argv) < 3:
    print(usage)
    sys.exit()
args = sys.argv[1:]
for arg in args:
    if len(arg.split('=')) == 1:
        print("Please connect arguments and values with '='")
        print(usage)
        sys.exit()
    option,value=arg.split('=')
    if option == '-B':   BEDGRAPH_IN  = value
    elif option == '-O': BEDGRAPH_OUT = value
    elif option == '-L': LENGTHS      = value
    elif option == '-S': STRAND       = value
    elif option == '-U': BED_UP       = value
    elif option == '-D': BED_DOWN     = value
if 'BEDGRAPH_IN' not in globals():
    print("ERROR: missing required argument -B=[input bedgraph]")
    print(usage)
    sys.exit()
if 'BEDGRAPH_OUT' not in globals():
    print("ERROR: missing required argument -O=[output bedgraph]")
    print(usage)
    sys.exit()
if 'LENGTHS' not in globals():
    print("ERROR: missing required argument -L=[lengths table]")
    print(usage)
    sys.exit()
if STRAND.lower() in ['+','plus','p']:
    STRAND='+'
elif STRAND.lower() in ['-','minus','m']:
    STRAND='-'
else:
    STRAND='.'
'''
'chromosomes' contains the lengths of all chromosomes the that BEDGRAPH contains values for.
Expects a two-column tab-separated file with:
    chromosome  length
Provided with the -L argument.
''' 
chromosomes={}
lengths_file=open(LENGTHS)
for line in lengths_file:
    chrom,length=line.rstrip().split('\t')
    chromosomes[chrom]=int(length)

mask_positions={}
for chrom in list(chromosomes.keys()):
    mask_positions[chrom]={}
    mask_positions[chrom]['-']=set()
    mask_positions[chrom]['+']=set()

if BED_UP != 'none':
    print('Upstream mask: '+BED_UP)
    mask_file=open(BED_UP)    
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

if BED_DOWN != 'none':
    print('Downstream mask: '+BED_DOWN)
    mask_file=open(BED_DOWN)    
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

bedgraph_file=open(BEDGRAPH_IN)
bedgraph_outfile=open(BEDGRAPH_OUT,'w')

for line in bedgraph_file:
    l=line.split('\t')
    chrom=l[0]
    end=int(l[1])
    if chrom not in chromosomes:
        continue
    if STRAND == 'both' and end not in mask_positions[chrom]['+'] and end not in mask_positions[chrom]['-']:
        bedgraph_outfile.write(line)
    elif end not in mask_positions[chrom][STRAND]:
        bedgraph_outfile.write(line)
bedgraph_outfile.close()
