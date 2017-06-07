import re,sys

usage="\
Outputs a column of count values representing the aggregate value in a BEDGRAPH file of all features in a GFF file\
\n\
Commandline arguments:\n\
    -B=[input bedgraph]           (required filepath)\n\
    -G=[input GFF]                (required filepath)\n\
    -L=[lengths table]            (required filepath; output from fasta_lengths.py)\n\
    -S=[strand]                   (default: both; options: both|. , plus|+ , minus|-)\n\
    -F=[feature type]             (default: gene; options: any string)\n\
    -V=[vertical normalization]   (default: norm; options: max,mean,norm)\n\
    -H=[horizontal normalization] (default: left; options: left,right,center,stretch<int>)\n\
\n\
This vector of values will grow in size over time based on the rule provided in -H:\n\
    left:         left-aligned\n\
    right:        right-aligned\n\
    center:       center-aligned\n\
    stretch<int>: All features will be resized into a vector of fixed length <int>\n\
\n\
Values added to the vector will be normalized prior according to -V:\n\
    max:          Sum of all values\n\
    mean:         Mean of all values\n\
    norm:         The values in each feature are scaled / max(abs(values))\n\
"

STRAND     = 'both'
FEATURE    = 'gene'
VERTICAL   = 'norm'
HORIZONTAL = 'left'

if len(sys.argv) < 2:
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
    elif option == '-G': GFF_IN      = value
    elif option == '-L': LENGTHS     = value
    elif option == '-S': STRAND      = value
    elif option == '-F': FEATURE     = value
    elif option == '-V': VERTICAL    = value
    elif option == '-H': HORIZONTAL  = value
if 'BEDGRAPH_IN' not in globals():
    print("ERROR: missing required argument -B=[input bedgraph]")
    print(usage)
    sys.exit()
if 'GFF_IN' not in globals():
    print("ERROR: missing required argument -G=[input GFF]")
    print(usage)
    sys.exit()
if 'LENGTHS' not in globals():
    print("ERROR: missing required argument -L=[lengths table]")
    print(usage)
    sys.exit()
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
for chrom,chromlen in chromosomes.items():
    coverage[chrom] = [0]*chromlen
coverage_file = open(BEDGRAPH_IN)
for line in coverage_file:
    chrom,start,end,count = line.rstrip().split()
    count = float(count)
    for i in range(int(start),int(end)):
        coverage[chrom][i] = count

'''
Begins reading the input GFF file and extracting the features of type -F.
'value_vector' is the set of running values that are aggregated.
The vector will grow in size over time based on the rule provided in -H:
    left:         left-aligned
    right:        right-aligned
    center:       center-aligned
    stretch<int>: All features will be resized into a vector of fixed length <int>
    
Values added to the vector will be normalized prior according to -V:
    max:          Sum of all values
    mean:         Mean of all values
    norm:         The values in each feature are scaled / max(abs(values))
'''
value_vector=[]
featurecount=0
gff=open(GFF_IN)
for line in gff:
    l=line.rstrip().split('\t')
    feature = l[2]
    if feature != FEATURE:
        continue
    chrom   = l[0]
    start   = int(l[3])-1
    end     = int(l[4])
    strand  = l[6]
    if start < 0 or end > chromosomes[chrom]:
        continue
    if strand == '+' and STRAND.lower() in ['-','minus','m']:
        continue
    if strand == '-' and STRAND.lower() in ['+','plus','p']:
        continue
    line_values = coverage[chrom][start:end]
    if strand == '-':
        line_values = list(reversed(line_values))
    if VERTICAL == 'norm':
        norm_value  = float(max([abs(i) for i in line_values]))
        if norm_value == 0:
            continue
        line_values = [float(i)/norm_value for i in line_values]
    veclen = len(value_vector)
    newlen = len(line_values)
    if 'stretch' in HORIZONTAL:
        print("Error: horizontal stretch not yet implemented.")
        sys.exit()
    elif HORIZONTAL == 'left':
        if veclen >= newlen:
            value_vector = [a+b for a,b in zip(value_vector,line_values+[0]*(veclen-newlen))]
        else:
            value_vector = [a+b for a,b in zip(value_vector+[0]*(newlen-veclen),line_values)]
    elif HORIZONTAL == 'right':
        if veclen >= newlen:
            value_vector = [a+b for a,b in zip(value_vector,[0]*(veclen-newlen)+line_values)]
        else:
            value_vector = [a+b for a,b in zip([0]*(newlen-veclen)+value_vector,line_values)]
    elif HORIZONTAL == 'center':
        print('Error: center-align not yet implemented')
        sys.exit()
    featurecount+=1
if VERTICAL in ['mean','norm']:
    value_vector = [i/featurecount for i in value_vector]

for v in value_vector:
    print('{:.10f}'.format(v))











