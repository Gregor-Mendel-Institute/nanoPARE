import sys,os,re,math
from collections import Counter

usage="\
Takes a BEDGRAPH of count values.\n\
All contiguous regions with values surpassing a given threshold\n\
are converted to features in an output BED file.\n\
\n\
Commandline arguments:\n\
    -B=[input bedgraph]           (required filepath)\n\
    -O=[output bed]               (required filepath)\n\
    -L=[lengths table]            (required filepath; output from fasta_lengths.py)\n\
    -T=[threshold]                (default: 0; number)\n\
    -M=[minimum feature length]   (default: 10; positive integer)\n\
    -V=[value type]               (default: sum; options: sum, max, mean)\n\
    -MV=[minimum value]           (default: 0; options: number)\n\
    -S=[strand]                   (default: both; options: both|. , plus|p|+ , minus|m|-)\n\
\n\
Any values immediately upstream of -U or downstream of -D will be set to zero.\n\
This can be used to mask sequence-specific artifacts, like TSO strand invasion or oligo-dT mispriming.\
"

THRESHOLD = 0
MINIMUM   = 10
VALUE     = 'sum'
MINVAL    = 0
STRAND    = 'both'

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
    if option == '-B':   BEDGRAPH_IN = value
    elif option == '-O': BED_OUT     = value
    elif option == '-L': LENGTHS     = value
    elif option == '-T': THRESHOLD   = float(value)
    elif option == '-M': MINIMUM     = int(value)
    elif option == '-V': VALUE       = value
    elif option == '-S': STRAND      = value
    elif option == '-MV': MINVAL     = float(value)
if 'BEDGRAPH_IN' not in globals():
    print("ERROR: missing required argument -B=[input bedgraph]")
    print(usage)
    sys.exit()
if 'BED_OUT' not in globals():
    print("ERROR: missing required argument -O=[output bedgraph]")
    print(usage)
    sys.exit()
if 'LENGTHS' not in globals():
    print("ERROR: missing required argument -L=[lengths table]")
    print(usage)
    sys.exit()
if type(THRESHOLD) not in [int,float]:
    print("ERROR: threshold must be numeric (-T=[threshold])")
    print(usage)
    sys.exit()
if type(MINIMUM) is not int:
    print("ERROR: minimum feature length must be an integer (-M=[minimum feature length])")
    print(usage)
    sys.exit()
VALUE=VALUE.lower()
if VALUE not in ['sum','max','mean']:
    print("ERROR: value type not supported (-V=[value type])")
    print(usage)
    sys.exit()
if STRAND.lower() in ['+','plus','p']:
    STRAND='+'
elif STRAND.lower() in ['-','minus','m']:
    STRAND='-'
else:
    STRAND='.'
    
def which(x,value=True):
    return [a for a,b in enumerate(x) if b==value]
    
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

'''
'coverage' is a dictionary of float vectors for each nucleotide in the genome.
Contains the value of the BEDGRAPH file at each position.
'''
coverage={}
thresh_coverage={}
for chrom,chromlen in chromosomes.items():
    coverage[chrom] = [0]*chromlen
    thresh_coverage[chrom] = [False]*chromlen

coverage_file = open(BEDGRAPH_IN)
for line in coverage_file:
    chrom,start,end,count = line.rstrip().split()
    count = float(count)
    for i in range(int(start),int(end)):
        coverage[chrom][i] = count
        thresh_coverage[chrom][i] = count > THRESHOLD

'''
Begins identifying continuous features, outputting a BED file.
Example line:
    Ath_chr1 1242364 1243615 thresh.358 0 + 31.1916 -1 -1 118
        chrom           chromosome name
        chromStart      leftmost end of the read (0-indexed)
        chromEnd        rightmost end of the read (0-indexed, open)
        name            unique name of each peak
        score           only used for visualizing data on a browser
        STRAND          directionality of peak (+ or -)
        signalValue     signal in peak, measured as VALUE type
        pValue          no pValue assigned, should be -1
        qValue          no qValue assigned, should be -1
        peak            position of point source for the peak (relative to chromStart, 0-indexed)
'''
featurecount = 0
outfile      = open(BED_OUT,'w')
for chrom,chromlen in sorted(list(chromosomes.items())):
    edges  = [int(thresh_coverage[chrom][0])]+[i-j for i,j in zip(thresh_coverage[chrom][1:],thresh_coverage[chrom][:chromlen])]
    starts = which(edges,1)
    ends   = which(edges,-1)
    if len(starts) > len(ends):
        ends += chromlen
    for chromStart,chromEnd in zip(starts,ends):
        if chromEnd-chromStart >= MINIMUM:
            featurecount += 1
            name = 'thresh.'+str(featurecount)
            coverage_subset = coverage[chrom][chromStart:chromEnd+1]
            if VALUE == 'sum':
                signalValue = sum(coverage_subset)
            elif VALUE == 'mean':
                signalValue = sum(coverage_subset)/len(coverage_subset)
            elif VALUE == 'max':
                signalValue = max(coverage_subset)
            if signalValue < MINVAL:
                continue
            if MINVAL > 0:
                signalValue = signalValue/MINVAL
            peak_positions = which(coverage_subset,max(coverage_subset))
            peak = int(sum(peak_positions)/len(peak_positions))
            outfile.write('\t'.join([str(i) for i in [chrom,chromStart,chromEnd,name,0,STRAND,signalValue,-1,-1,peak]])+'\n')
outfile.close()
