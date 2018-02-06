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

value_vector=[]

featurecount=0

if 'stretch' in args.HORIZONTAL:
    args.HORIZONTAL = int(args.HORIZONTAL.replace('stretch',''))
    value_vector = [float(0)]*args.HORIZONTAL

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
    
    if args.FEATURE != 'transcript':
        sequence = ''.join([genome[chrom][i] for i in positions])
        if strand == '-':
            aa,ss,f = fu.longest_orf(fu.rc(sequence))
            ORFstart,ORFstop = ss
            if args.FEATURE == '5UTR':
                positions = positions[-ORFstart:]
            elif args.FEATURE == 'CDS':
                positions = positions[-ORFstop:-ORFstart]
            elif args.FEATURE == '3UTR':
                positions = positions[:-ORFstop]
        else:
            aa,ss,f = fu.longest_orf(sequence)
            ORFstart,ORFstop = ss
            if args.FEATURE == '5UTR':
                positions = positions[:ORFstart]
            elif args.FEATURE == 'CDS':
                positions = positions[ORFstart:ORFstop]
            elif args.FEATURE == '3UTR':
                positions = positions[ORFstop:]
    
    
    if len(positions) == 0:
        continue
    
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
    

    
    if args.ANTISENSE:
        if strand == '+':
            line_values = values_minus
        elif strand == '-':
            line_values = values_plus
    else:
        if strand == '+':
            line_values = values_plus
        elif strand == '-':
            line_values = values_minus
    if args.NONSTRANDED:
        line_values = [a+b for a,b in zip(values_ns,line_values)]
    
    if not line_values:
        continue
    
    if args.VERTICAL == 'norm':
        norm_value  = float(max([abs(i) for i in line_values]))
        if norm_value == 0:
            continue
        line_values = [float(i)/norm_value for i in line_values]
    
    if type(args.HORIZONTAL) is int:
        bins = [float(0)]*args.HORIZONTAL
        partition = float(len(line_values))/args.HORIZONTAL
        for i in range(0,args.HORIZONTAL):
            values_to_bin = line_values[int(round(partition*(i))):int(round(partition*(i+1)))]
            if not len(values_to_bin)==0:
                if args.VERTICAL == 'max':
                    bins[i] = sum(values_to_bin)
                else:
                    bins[i] = sum(values_to_bin)/len(values_to_bin)
        value_vector = [a+b for a,b in zip(value_vector,bins)]
    else:
        veclen = len(value_vector)
        newlen = len(line_values)
        if args.HORIZONTAL == 'left':
            if veclen >= newlen:
                value_vector = [a+b for a,b in zip(value_vector,line_values+[0]*(veclen-newlen))][:args.MAXLEN]
            else:
                value_vector = [a+b for a,b in zip(value_vector+[0]*(newlen-veclen),line_values)][:args.MAXLEN]
        elif args.HORIZONTAL == 'right':
            if veclen >= newlen:
                value_vector = [a+b for a,b in zip(value_vector,[0]*(veclen-newlen)+line_values)][-args.MAXLEN:]
            else:
                value_vector = [a+b for a,b in zip([0]*(newlen-veclen)+value_vector,line_values)][-args.MAXLEN:]
        elif args.HORIZONTAL == 'center':
            print('Error: center-align not yet implemented')
            sys.exit()
    featurecount+=1
    if args.SIGNAL:
        if not signal_out:
            signal_out = open(args.SIGNAL,'w')
        if sum(line_values) > 0:
            signal_out.write('{}\n'.format(ID))

if args.SIGNAL:
    signal_out.close()

if args.VERTICAL in ['mean','norm']:
    value_vector = [i/featurecount for i in value_vector]

if args.VERTICAL.lower() == 'rpm':
    value_vector = [i/million_counts for i in value_vector]

for v in value_vector:
    print('{:.10f}'.format(v))








