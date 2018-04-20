import argparse
import fasta_utils as fu

###################
# INPUT ARGUMENTS #
###################

desc = (
    "Takes a BED file of features and finds read coverage."
    "Imports read coverage for all reads and for those with"
    "upstream untemplated G (uuG)."
    "Separates the BED file into high-uuG (capped) and low-uuG (uncapped)."
    "\n"
)

parser = argparse.ArgumentParser(description=desc)

# add arguments to the ArgumentParser
parser.add_argument(
    'bed_in', metavar='FEATURES_BED',
    type=str, help='input bed file', action='store'
)
parser.add_argument(
    'all_bedgraphs', metavar='ALL_BEDGRAPH',
    type=str, help='input bedgraph files (+/-) of all end reads',
    action='store', nargs=2
)
parser.add_argument(
    'uug_bedgraphs', metavar='UUG_BEDGRAPH',
    type=str, help='input bedgraph files (+/-) of uuG-containing reads',
    action='store', nargs=2
)
parser.add_argument(
    '-C', '--capped_output', dest='capped_output', metavar='"file"',
    type=str, help='output bed for capped features', action='store',
    default='capped_features.bed'
)
parser.add_argument(
    '-U', '--uncapped_output', dest='uncapped_output', metavar='"file"',
    type=str, help='output bed for uncapped features', action='store',
    default='uncapped_features.bed'
)
parser.add_argument(
    '--minimum', dest='minimum', type=int,
    help='minumum number of reads in a feature to keep it',
    default=1
)
parser.add_argument(
    '--min_uug', dest='min_uug', type=float,
    help='minumum %uuG to consider capped',
    default=None
)


args = parser.parse_args()

####################
# DEFINE FUNCTIONS #
####################

def which(x,value=True):
    """Returns a list of locations in x that satisfy value"""
    return [a for a,b in enumerate(x) if b==value]


def table(x):
    y = sorted(list(set(x)))
    z = {}
    for i in y:
        z[i] = len([1 for b in x if b == i])
    return z


def percentile_bin(values, normalize=False, digits=2):
    """Returns an ordered matrix of bins that represent a histogram
    of values from 0-1 by an increment of 10^(-digits).
    Each bin contains the number of elements 'v' in 'values' that satisfy:
    
    [bin - 1] <= v < [bin]
    """
    outputbin = list(range(0,10**digits + 1))
    
    if normalize:
        minval = min(values)
        maxval = max(values) - minval
        values = [(float(i) - minval)/maxval for i in values]
    
    valuetable = table(values)
    allvalues = list(valuetable.keys())
    
    minval = 0
    maxval = 1
    for b in range(len(outputbin)):
        lowval = float(b)/(10**digits)
        highval = float(b + 1)/(10**digits)
        outputbin[b] = sum(
                [
                    valuetable[i]
                    for i in [
                        j
                        for j in allvalues
                        if j >= lowval and j < highval
                    ]
                ]
        )
    
    return outputbin


####################
# LOAD ENVIRONMENT #
####################

# 'features' contains all elements of the input BED file
features = {}
for strand in ['+', '-']:
    features[strand] = {}

feature_file = open(args.bed_in)
for line in feature_file:
    if line[0] == '#':
        continue
    
    l = line.rstrip().split('\t')
    score = l[4]
    chrom = l[0]
    start_pos = int(l[1])
    end_pos = int(l[2])
    # peak_pos = int(l[6]) + start_pos
    if len(l) > 6:
        secondary = l[6:]
    else:
        secondary = None
    
    strand = l[5]
    readname = l[3]

    if chrom not in features[strand]:
        features[strand][chrom] = {}
    
    features[strand][chrom][start_pos] = [
        (start_pos,end_pos),
        score,
        readname,
        secondary
    ]

feature_file.close()

# 'coverage' contains dictionaries
# of float vectors for each nucleotide in the genome.
# Saves the value of the respective BEDGRAPH file at each position.
coverage = {}
coverage['all'] = {}
coverage['all']['+'] = {}
coverage['all']['-'] = {}
coverage['uug'] = {}
coverage['uug']['+'] = {}
coverage['uug']['-'] = {}
    
all_plus = args.all_bedgraphs[0]
all_minus = args.all_bedgraphs[1]
uug_plus = args.uug_bedgraphs[0]
uug_minus = args.uug_bedgraphs[1]

print('Loading all plus...')
for line in open(all_plus):
    if line[0] =='#':
        continue
    
    chrom,start,end,count = line.rstrip().split()
    count = float(count)
    if chrom not in coverage['all']['+']:
        coverage['all']['+'][chrom] = {}
    
    for pos in range(int(start),int(end)):
        coverage['all']['+'][chrom][pos] = count

print('Loading all minus...')
for line in open(all_minus):
    if line[0] =='#':
        continue
    
    chrom,start,end,count = line.rstrip().split()
    if chrom not in coverage['all']['-']:
        coverage['all']['-'][chrom] = {}
    
    count = float(count)
    for pos in range(int(start),int(end)):
        coverage['all']['-'][chrom][pos] = count

print('Loading uuG plus...')
for line in open(uug_plus):
    if line[0] =='#':
        continue
    
    chrom,start,end,count = line.rstrip().split()
    if chrom not in coverage['uug']['+']:
        coverage['uug']['+'][chrom] = {}
    
    count = float(count)
    for pos in range(int(start),int(end)):
        coverage['uug']['+'][chrom][pos] = count

print('Loading uuG minus...')
for line in open(uug_minus):
    if line[0] =='#':
        continue
    
    chrom,start,end,count = line.rstrip().split()
    if chrom not in coverage['uug']['-']:
        coverage['uug']['-'][chrom] = {}
    
    count = float(count)
    for pos in range(int(start),int(end)):
        coverage['uug']['-'][chrom][pos] = count

print("Evaluating features...")
allfeatures = []
chromosomes = set()

for strand in ['+','-']:
    for key in features[strand].keys():
        chromosomes.add(key)


for strand in ['+','-']:
    for chrom in sorted(list(chromosomes)):
        if chrom in features[strand]:
            for currentpos,feature in features[strand][chrom].items():
                # read feature traits from the dict
                startend,score,readname,secondary = feature
                start,end = [int(i) for i in startend]
                
                # Get coverage values overlapping this feature
                all_reads = float(0)
                if chrom in coverage['all'][strand]:
                    for pos in range(start,end):
                        all_reads += coverage['all'][strand][chrom].get(pos,0)
                
                uug_reads = 0
                if chrom in coverage['uug'][strand]:
                    for pos in range(start,end):
                        uug_reads += coverage['uug'][strand][chrom].get(pos,0)
                
                # discard a feature with less than minimum reads
                if all_reads < args.minimum:
                    continue

                # Calculate percent uuG and pipe the feature appropriately
                percent_uuG = float(uug_reads)/all_reads
                
                # Make a BED feature string
                feature_to_write = '\t'.join(
                    [
                        chrom,
                        str(start),
                        str(end),
                        readname,
                        str(score),
                        strand,
                        '\t'.join(secondary),
                        str(round(percent_uuG,3))
                    ]
                ) + '\n'
                
                allfeatures.append((feature_to_write,percent_uuG))

if args.min_uug:
    cutoff = args.min_uug
    print("# Cutoff: {}".format(cutoff))
else:
    # Determine the cutoff for partitioning capped and noncapped features
    # by locating the local minimum in distribution of uuG reads
    print("Determining percent uuG cutoff...")
    all_uuG = [b for a,b in allfeatures]
    uug_dist = percentile_bin(all_uuG, normalize=False, digits=2)
    uug_sub = uug_dist[:int(len(uug_dist)/4)]
    minima = which(uug_sub,min(uug_sub))
    if len(minima) == 1:
        cutoff = float(minima[0])/100
    else:
        print("# WARNING: More than 1 local minimum")
        cutoff = float(sum(minima)/len(minima))/100
    print("# Local minimum: {}".format(cutoff))
    print("# uuG distribution:")
    for a,b in zip(range(0,10**2 + 1),uug_dist):
        print('{}\t{}'.format(a,b))


# Write the features to 'capped' or 'noncapped' output based on the cutoff
capped_out = open(args.capped_output, 'w')
uncapped_out = open(args.uncapped_output, 'w')
for feature_to_write,percent_uuG in allfeatures:
    if percent_uuG < cutoff:
        uncapped_out.write(feature_to_write)
    else:
        capped_out.write(feature_to_write)

capped_out.close()
uncapped_out.close()
