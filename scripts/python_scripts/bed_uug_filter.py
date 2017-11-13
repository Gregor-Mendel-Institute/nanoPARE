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
    'genome', metavar='FASTA',
    type=str, help='genome FASTA file'
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
    default=5
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
    outputbin = range(0,10**digits + 1)
    
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

# 'genome' is a dict of strings of the genome: {name:string}
genome = fu.import_genome(args.genome)

# 'chromosomes' contains the lengths of all chromosomes
chromosomes = {}
for chrom in genome.keys():
    length = len(genome[chrom])
    chromosomes[chrom] = length

# 'features' contains all elements of the input BED file
features = {}
for readtype in ['TSS', 'TES']:
    for strand in ['plus', 'minus']:
        features[readtype + strand] = {}

for chromosome in chromosomes.keys():
    for readtype in ['TSS', 'TES']:
        for strand in ['plus', 'minus']:
            features[readtype + strand][chromosome] = {}

feature_file = open(args.bed_in)
for line in feature_file:
    if line[0] == '#':continue
    l = line.rstrip().split('\t')
    score = l[4]
    chrom = l[0]
    start_pos = int(l[1])
    end_pos = int(l[2])
    peak_pos = int(l[6]) + start_pos
    secondary = l[-1]
    
    readname = l[3]
    readtype = 'TSS'
    if l[5] == '+':
        strand = 'plus'
    elif l[5] == '-':
        strand = 'minus'
    
    features[readtype + strand][chrom][peak_pos] = [
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
coverage['all']['plus'] = {}
coverage['all']['minus'] = {}
coverage['uug'] = {}
coverage['uug']['plus'] = {}
coverage['uug']['minus'] = {}
for chrom,chromlen in chromosomes.items():
    coverage['all']['plus'][chrom] = [0]*chromlen
    coverage['all']['minus'][chrom] = [0]*chromlen
    coverage['uug']['plus'][chrom] = [0]*chromlen
    coverage['uug']['minus'][chrom] = [0]*chromlen    
    
all_plus = args.all_bedgraphs[0]
all_minus = args.all_bedgraphs[1]
uug_plus = args.uug_bedgraphs[0]
uug_minus = args.uug_bedgraphs[1]

print('Loading all plus...')
for line in open(all_plus):
    chrom,start,end,count = line.rstrip().split()
    count = float(count)
    coverage['all']['plus'][chrom][int(start):int(end)] = [count]*(int(end)-int(start))

print('Loading all minus...')
for line in open(all_minus):
    chrom,start,end,count = line.rstrip().split()
    count = float(count)
    coverage['all']['minus'][chrom][int(start):int(end)] = [count]*(int(end)-int(start))

print('Loading uuG plus...')
for line in open(uug_plus):
    chrom,start,end,count = line.rstrip().split()
    count = float(count)
    coverage['uug']['plus'][chrom][int(start):int(end)] = [count]*(int(end)-int(start))

print('Loading uuG minus...')
for line in open(uug_minus):
    chrom,start,end,count = line.rstrip().split()
    count = float(count)
    coverage['uug']['minus'][chrom][int(start):int(end)] = [count]*(int(end)-int(start))

print("Evaluating features...")
allfeatures = []
for readtype in ['TSS']:
    for strand in ['plus','minus']:
        for chrom in sorted(list(chromosomes.keys())):
            for peakpos,feature in features[readtype + strand][chrom].items():
                # read feature traits from the dict
                startend,score,readname,secondary = feature
                start,end = [int(i) for i in startend]
                
                # slice coverage dicts for values overlapping this feature
                all_reads = coverage['all'][strand][chrom][start:end]
                uug_reads = coverage['uug'][strand][chrom][start:end]
                
                # get the nucleotide sequence immediately upstream
                if strand == 'plus':
                    upstream_seqs = genome[chrom][(start-1):(end-1)]
                elif strand == 'minus':
                    upstream_seqs = reversed(
                        fu.rc(
                            genome[chrom][(start+1):(end+1)]
                        )
                    )
                
                # discard a feature with less than minimum reads
                if sum(all_reads) < args.minimum:
                    continue
                
                # get positions without a templated upstream G
                non_tuG = which([i != 'G' for i in upstream_seqs])
                
                # calculate the total reads and uuG reads at non-tuG sites
                informative_all = sum(
                    [
                        all_reads[i]
                        for i in non_tuG
                    ]
                )
                
                # Discard peaks with no reads informative of
                # whether or not the feature is capped
                if informative_all == 0:
                    continue
                
                informative_uuG = sum(
                    [
                        uug_reads[i]
                        for i in non_tuG
                    ]
                )
                
                # Calculate percent uuG and pipe the feature appropriately
                percent_uuG = float(informative_uuG)/informative_all
                
                # Make a BED feature string
                if strand == 'plus':
                    s = '+'
                elif strand == 'minus':
                    s = '-'
                
                feature_to_write = '\t'.join(
                    [
                        chrom,
                        str(start),
                        str(end),
                        readname,
                        str(score),
                        s,
                        str(int(peakpos)-start),
                        secondary,
                        str(round(percent_uuG,3))
                    ]
                ) + '\n'
                
                allfeatures.append((feature_to_write,percent_uuG))


# Determine the cutoff for partitioning capped and noncapped features
# by locating the local minimum in distribution of uuG reads
print("Determining percent uuG cutoff...")
all_uuG = [b for a,b in allfeatures]
uug_dist = percentile_bin(all_uuG, normalize=False, digits=2)
uug_sub = uug_dist[:len(uug_dist)/4]
minima = which(uug_sub,min(uug_sub))
if len(minima) == 1:
    cutoff = float(minima[0])/100
else:
    print("# WARNING: More than 1 local minimum")
    cutoff = float(sum(minima)/len(minima))/100

print("# Cutoff: {}".format(cutoff))
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
