import sys,os,re,math
from collections import Counter

usage="\
Takes a BEDGRAPH of count values.\n\
Each point is converted to a density distribution and all densities are summed\n\
for the output model.\n\
\n\
Commandline arguments:\n\
    -B=[input bedgraph]           (required filepath)\n\
    -O=[output bedgraph]          (required filepath)\n\
    -L=[lengths table]            (required filepath; output from fasta_lengths.py)\n\
    -K=[kernel function]          (default: gaussian; options: gaussian|normal, laplace|folded, uniform , triangle)\n\
    -H=[bandwidth]                (default: 50; options: integer)\n\
    -S=[sigma value]              (default: 5; options: 1-6)\n\
    -D=[significant digits]       (default: 4; options: integer)\n\
    -P=[only positive values]     (default: FALSE)\n\
\n\
Any values immediately upstream of -U or downstream of -D will be set to zero.\n\
This can be used to mask sequence-specific artifacts, like TSO strand invasion or oligo-dT mispriming.\
"

KERNEL    = 'laplace'
BANDWIDTH = 50
SIGMA     = 5
DIGITS    = 3
ONLY_POSITIVE = False

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
    if option == '-B':   BEDGRAPH_IN   = value
    elif option == '-O': BEDGRAPH_OUT  = value
    elif option == '-L': LENGTHS       = value
    elif option == '-K': KERNEL        = value
    elif option == '-H': BANDWIDTH     = int(value)
    elif option == '-S': SIGMA         = int(value)
    elif option == '-D': DIGITS        = int(value)
    elif option == '-P': ONLY_POSITIVE = bool(value)
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
if SIGMA not in [1,2,3,4,5,6]:
    print("ERROR: sigma must be an integer 1-6 (-S=[sigma value])")
    print(usage)
    sys.exit()
if type(BANDWIDTH) is not int:
    print("ERROR: bandwidth must be an integer (-H=[bandwidth])")
    print(usage)
    sys.exit()
if type(DIGITS) is not int or DIGITS <= 0:
    print("ERROR: significant must be a positive integer (-D=[significant digits])")
    print(usage)
    sys.exit()
'''
Gives a density distribution as a list of values centered symmetrically around 'center',
stretching 'bandwidth' units left and right and multiplying the density by 'value' (e.g. raw or normalized counts).
The density function to use is defined by 'd_fun'. Options:
    Gaussian - 'bandwidth' is 1-sigma in a Gaussian distribution function; returns a list of values +- 3 sigma in length.
    Uniform  - Every position in [center-bandwidth : center+bandwidth+1] is given an equal value such that sum(positions) = 1
    Triangle - Center-weighted symmetric triangular distribution.
'''

def gauss(x,s):
    return math.exp(-float(x)**2/(2.0*float(s)**2))/(math.sqrt(2.0*math.pi)*s)
    
def laplace(x,s,m=0):
    return math.exp(-abs(float(x)-m)/s)/(2*s)
    
def which(x,value=True):
    return [a for a,b in enumerate(x) if b==value]
    
def notwhich(x,value=0):
    return [a for a,b in enumerate(x) if b!=value]
    
if KERNEL.lower() in ['gauss','gaussian','normal','norm']:
    k = [gauss(i,BANDWIDTH) for i in range(-BANDWIDTH*SIGMA,BANDWIDTH*SIGMA+1)]
    k_len = len(k)
elif KERNEL.lower() in ['laplace','fold','folded','mountain']:
    k = [laplace(i,BANDWIDTH) for i in range(-BANDWIDTH*SIGMA,BANDWIDTH*SIGMA+1)]
    k_len = len(k)
elif KERNEL.lower() in ['uniform','flat']:
    k_len = len(range(-BANDWIDTH,BANDWIDTH+1))
    k = [1.0/k_len]*k_len
elif KERNEL.lower() in ['triangle','triangular']:
    k_len = len(range(-BANDWIDTH,BANDWIDTH+1))
    k = [BANDWIDTH-abs(BANDWIDTH-i) for i in range(k_len)]
    k_sum = sum(k)
    k = [float(i)/k_sum for i in k]
else:
    print("ERROR: kernel argument not supported (-K=[kernel function])")
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
'smoothed_coverage' is an empty dictionary with the same structure as 'coverage',
to be filled with smoothed values from the kernel function.
'''
coverage={}
smoothed_coverage={}
for chrom,chromlen in chromosomes.items():
    coverage[chrom] = [0]*chromlen
    smoothed_coverage[chrom] = [0]*chromlen

coverage_file = open(BEDGRAPH_IN)
for line in coverage_file:
    chrom,start,end,count = line.rstrip().split()
    count = float(count)
    for i in range(int(start),int(end)):
        coverage[chrom][i] = count

'''
Execute kernel smoothing:
    For every chromosome listed in the lengths table,
    transform each nonempty value in 'coverage' to a distribution centered at its position.
    The distribution is defined by the kernel function above, multiplied by the value of 'coverage'
    Add the values of the distribution to the respective positions in 'smoothed_coverage'
'''
outfile=open(BEDGRAPH_OUT,'w')
for chrom,chromlen in sorted(list(chromosomes.items())):
    nonempty=[(loc,val) for loc,val in zip(range(chromlen),coverage[chrom]) if val]
    for loc,val in nonempty:
        k_locs = range(loc-(k_len-1)/2,loc+(k_len-1)/2+1)
        k_vals = [i*val for i in k]
        for l,v in zip(k_locs,k_vals):
            if l<0 or l>=chromlen:
                continue
            smoothed_coverage[chrom][l] += v
    del coverage[chrom]
    if ONLY_POSITIVE:
        smooth_round = [round(i,DIGITS) if round(i,DIGITS)>0 else False for i in smoothed_coverage[chrom]]
    else:
        smooth_round = [round(i,DIGITS) if round(i,DIGITS)!=0 else False for i in smoothed_coverage[chrom]]
    del smoothed_coverage[chrom]
    edges = [smooth_round[0]]+[i-j for i,j in zip(smooth_round[1:],smooth_round[:chromlen])]
    diff  = notwhich(edges)
    del edges
    diffpairs = zip(diff,diff[1:])
    del diff
    for locs,val in zip(diffpairs,[smooth_round[i[0]] for i in diffpairs]):
        if val != 0:
            outfile.write('\t'.join([chrom,str(locs[0]),str(locs[1]),str(val)])+'\n')
outfile.close()
