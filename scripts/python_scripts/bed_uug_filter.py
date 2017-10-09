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
    'bed_in', metavar='BED',
    type=str, help='input bed file', action='store'
)
parser.add_argument(
    'all_bedgraphs', metavar='BEDGRAPH',
    type=str, help='input bedgraph files (+/-) of all end reads',
    action='store', nargs=2
)
parser.add_argument(
    'uug_bedgraphs', metavar='BEDGRAPH',
    type=str, help='input bedgraph files (+/-) of uuG-containing reads',
    action='store', nargs=2
)
parser.add_argument(
    'genome', metavar='FASTA',
    type=str, help='table of chromosome lengths'
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
parser.add_argument(
    '--cutoff', dest='cutoff', type=float,
    help='lowest proportion uuG to be considered capped',
    default=0.05
)

args = parser.parse_args()

####################
# DEFINE FUNCTIONS #
####################

def which(x,value=True):
    """Returns a list of locations in x that satisfy value"""
    return [a for a,b in enumerate(x) if b==value]


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
    l = line.rstrip().split()
    score = float(l[4])
    chrom = l[0]
    start_pos = int(l[1])
    end_pos = int(l[2])
    peak_pos = int(l[6]) + start_pos
    if len(l) == 8:
        if l[7] == '-':
            secondary_peaks = []
        else:
            secondary_peaks = [
                int(i) + start_pos
                for i in l[7].split(',')
            ]
    else:
        secondary_peaks = []
    
    readname = l[3]
    readtype = l[3].split('.')[0]
    strand = l[3].split('.')[1]
    
    features[readtype + strand][chrom][peak_pos] = [
        (start_pos,end_pos),
        score,
        readname,
        secondary_peaks
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

capped_out = open(args.capped_output, 'w')
uncapped_out = open(args.uncapped_output, 'w')

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
                        ','.join([str(i - start) for i in secondary]),
                        str(round(percent_uuG,3))
                    ]
                ) + '\n'
                
                if percent_uuG < args.cutoff:
                    uncapped_out.write(feature_to_write)
                else:
                    capped_out.write(feature_to_write)

capped_out.close()
uncapped_out.close()
