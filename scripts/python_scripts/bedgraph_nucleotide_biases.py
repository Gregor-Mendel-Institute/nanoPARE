import re
import sys
import argparse
import gff_utils as gu
import fasta_utils as fu
import bedgraph_utils as bu

######################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-A','--annotation',dest='ANNOTATION',
                    help="Path to a reference annotation GTF/GFF.",
                    required=True)
parser.add_argument('-G','--genome',dest='GENOME',
                    help="Path to a reference genome FASTA.")
parser.add_argument("--subset", dest='SUBSET',
                    help="Subset of feature IDs to use from the annotation.", 
                    default=None, type=str)
parser.add_argument('-P',dest='PLUS_BEDGRAPH',
                    help="Path to nonstranded bedgraph file(s).",
                    default=[], nargs = '+')
parser.add_argument('-M',dest='MINUS_BEDGRAPH',
                    help="Path to nonstranded bedgraph file(s).",
                    default=[], nargs = '+')
parser.add_argument('-N',dest='NONSTRANDED_BEDGRAPH',
                    help="Path to nonstranded bedgraph file(s).",
                    default=[], nargs = '+')
parser.add_argument('-F','--feature',dest='FEATURE',
                    help="Feature type to plot",
                    default='transcript',choices=['CDS','transcript','5UTR','3UTR'],type=str)
parser.add_argument('-V','--vertical',dest='VERTICAL',
                    help="Vertical normalization method",
                    default='norm',type=str)
parser.add_argument('-H','--horizontal',dest='HORIZONTAL',
                    help="Horizontal normalization method",
                    default='stretch100',type=str)
parser.add_argument('-B','--buffer', dest='BUFFER',
                    help="Number of nucleotides to add to each end",
                    default=0, type=int)
parser.add_argument('--maxlen', dest='MAXLEN',
                    help="Longest the vector can grow.",
                    default=0, type=int)
parser.add_argument("--features_with_signal", dest='SIGNAL',
                    help="Filepath to output gene IDs with signal.", 
                    default=None, type=str)
parser.add_argument("--offset", dest='OFFSET',
                    help="Number of nucleotides up/down to shift the pattern search", 
                    default=1, type=int)
parser.add_argument("--length", dest='LENGTH',
                    help="length of sequence motif to analyze", 
                    default=4, type=int)
parser.add_argument("--splice", dest='SPLICE',
                    help="Removes intronic sequences.", 
                    default=False, action='store_true')
parser.add_argument("--antisense", dest='ANTISENSE',
                    help="Reports metaplot antisense to the annotation.", 
                    default=False, action='store_true')
parser.add_argument("--ns", dest='NONSTRANDED',
                    help="Reports nonstranded data.", 
                    default=False, action='store_true')
parser.add_argument("--no_doubles", dest='NO_DOUBLES',
                    help="Don't count any nucleotide more than once.", 
                    default=False, action='store_true')
args = parser.parse_args()
######################################################################

def which(x,value=True):
    """Returns a list of locations in x that satisfy value"""
    return [a for a,b in enumerate(x) if b==value]


def notwhich(x,value=0):
    """Returns a list of locations in x that do not satisty value"""
    return [a for a,b in enumerate(x) if b!=value]


def flatten(list_of_lists):
    """Collapses a list/tuple of lists into a single list"""
    return [item for sublist in list_of_lists for item in sublist]


if args.GENOME:
    genome = fu.import_genome(args.GENOME)
else:
    if args.FEATURE != 'transcript':
        print("ERROR: cannot locate {} features without a reference genome.".format(args.FEATURE))
        print("Provide genome FASTA file with -G")
        sys.exit(1)

coverage = None
if args.PLUS_BEDGRAPH:
    for i in args.PLUS_BEDGRAPH:
        if not coverage:
            coverage = bu.parse_bedgraph(i,'+')
        else:
            bu.add_bedgraph(coverage,i,'+')
    
if args.MINUS_BEDGRAPH:
    for i in args.MINUS_BEDGRAPH:
        if not coverage:
            coverage = bu.parse_bedgraph(i,'-')
        else:
            bu.add_bedgraph(coverage,i,'-')

if args.NONSTRANDED_BEDGRAPH:
    for i in args.NONSTRANDED_BEDGRAPH:
        if not coverage:
            coverage = bu.parse_bedgraph(i,'.')
        else:
            bu.add_bedgraph(coverage,i,'.')

if not coverage:
    print("ERROR: No bedgraph files parsed.")
    sys.exit(1)

if args.VERTICAL.lower() == 'rpm':
    million_counts = sum(
        [
            sum(coverage['+'][i].values())
            for i in coverage['+'].keys()
        ] + [
            sum(coverage['-'][i].values())
            for i in coverage['-'].keys()
        ] + [
            sum(coverage['.'][i].values())
            for i in coverage['.'].keys()
        ]
    )*(10**-6)

motifs_total = {}
motifs_positive = {}

ref_filetype = args.ANNOTATION.split('.')[-1].lower()
if ref_filetype == 'gff':
    ref_filetype = 'gff3'

ref_transcripts = gu.get_file_content(args.ANNOTATION,ref_filetype)
ref_IDs = sorted(list(ref_transcripts.keys()))
print('# {} reference features: {}'.format(
    len(ref_IDs),
    args.ANNOTATION
))

if args.SUBSET:
    sub_IDs = [l.rstrip() for l in open(args.SUBSET)]
    ref_set = set(ref_IDs)
    picked_IDs = [i for i in sub_IDs if i in ref_set]
    print('# {} picked features: {}'.format(
        len(picked_IDs),
        args.SUBSET
    ))
else:
    picked_IDs = ref_IDs

signal_out = None
featurecount = 0

for ID in picked_IDs:
    # Get feature attributes
    feature = ref_transcripts[ID]
    chrom = feature.chrom
    strand = feature.strand
    start = int(feature.start)
    end = int(feature.stop)
    start = start - 1
    end = end
    
    if args.SPLICE:
        exon_starts = feature.get_exon_start()
        exon_ends = feature.get_exon_end()
        if len(exon_starts) != len(exon_ends):
            print('# WARNING: inconsistent exon structure with {}'.format(ID))
            continue
        
        positions = flatten(
            [list(range(a - 1,b)) for a,b in zip(exon_starts,exon_ends)]
        )
    else:
        positions = list(range(start,end))
    
    if strand == '-':
        positions = positions[::-1]
    
    length = len(positions)
    
    sequence = ''.join([genome[chrom][i] for i in positions])
    if strand == '-':
        sequence = fu.complement(sequence)
    
    aa,ss,f = fu.longest_orf(sequence)
    
    ORFstart,ORFstop = ss
    if args.FEATURE == '5UTR':
        positions = positions[:ORFstart]
        sequence = sequence[:ORFstart]
    elif args.FEATURE == 'CDS':
        positions = positions[ORFstart:ORFstop]
        sequence = sequence[ORFstart:ORFstop]
    elif args.FEATURE == '3UTR':
        positions = positions[ORFstop:]
        sequence = sequence[ORFstop:]
    
    if args.BUFFER:
        if strand == '+':
            positions = list(range(positions[0] - args.BUFFER, positions[0])) + positions + list(range(positions[-1] + 1, positions[-1] + 1 + args.BUFFER))
        else:
            positions = list(range(positions[0] + args.BUFFER, positions[0], -1)) + positions + list(range(positions[-1] - 1, positions[-1] - 1 - args.BUFFER, -1))
    
    if start < 0:
        continue
    
    if chrom in coverage['+']:
        values_plus = [coverage['+'][chrom].get(i,0) for i in positions]
        if args.NO_DOUBLES:
            for i in positions:
                if i in coverage['+'][chrom]:
                    coverage['+'][chrom][i] = 0
    else:
        values_plus = [0 for i in positions]
    
    if chrom in coverage['-']:
        values_minus = [coverage['-'][chrom].get(i,0) for i in positions]
        if args.NO_DOUBLES:
            for i in positions:
                if i in coverage['-'][chrom]:
                    coverage['-'][chrom][i] = 0
    else:
        values_minus = [0 for i in positions]
    
    if chrom in coverage['.']:
        values_ns = [coverage['.'][chrom].get(i,0) for i in positions]
        if args.NO_DOUBLES:
            for i in positions:
                if i in coverage['.'][chrom]:
                    coverage['.'][chrom][i] = 0
    else:
        values_ns = [0 for i in positions]
    
    if args.NONSTRANDED:
        line_values = values_ns
    elif args.ANTISENSE:
        sequence = fu.complement(sequence)
        if strand == '+':
            line_values = values_minus
        elif strand == '-':
            line_values = values_plus
    else:
        if strand == '+':
            line_values = values_plus
        elif strand == '-':
            line_values = values_minus
        
    # Iterate through all positions of the transcript and calculate signal coverage
    for i in range(len(sequence)):
        subseq = sequence[(i+args.OFFSET):(i+args.OFFSET+args.LENGTH)]
        if len(subseq) != args.LENGTH:
            continue
        motifs_total[subseq] = motifs_total.get(subseq,0) + 1
        signal = line_values[i]
        if signal > 0:
            motifs_positive[subseq] = motifs_positive.get(subseq,0) + 1
    
    featurecount+=1

positions_examined = sum(motifs_total.values())
positions_with_signal = sum(motifs_positive.values())
for motif in sorted(motifs_total.keys()):
    motif_proportion = float(motifs_total[motif]) / positions_examined
    positive_proportion = float(motifs_positive.get(motif,0)) / positions_with_signal
    print('{}\t{}'.format(
        motif,
        positive_proportion/motif_proportion
        )
    )








