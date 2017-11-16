import re
import sys
import argparse
import gff_utils as gu

desc = (
    "Calculates a scaling ratio between two bedgraph files."
    " Each gene in a GTF/GFF file is examined, and a mean read depth"
    " is calculated for every position where reads are found in the gene"
    " for both a 'positive' and 'negative' bedgraph file."
)

parser = argparse.ArgumentParser(description=desc)

# add arguments to the ArgumentParser
parser.add_argument(
    '-P', '--pos_bg', dest='pos_bg', type=str, 
    help='input positive bedgraph filepath', required=True
)
parser.add_argument(
    '-N', '--neg_bg', dest='neg_bg', type=str, 
    help='input negative bedgraph filepath', required=True
)
parser.add_argument(
    '-A', '--anno', dest='reference_GFF', type=str, 
    help='input annotation GFF/GTF filepath', required=True
)
parser.add_argument(
    '-S', '--strand', dest='strand', type=str, 
    help='(+/-) which strand coverage data is on', required=True
)
parser.add_argument(
    '--flanking', dest='flanking', type=int, 
    help='number of nucleotides to extend annotations', default=200
)
parser.add_argument(
    '--maxband', dest='maxband', type=int, 
    help='maximum bandwidth', default=50
)
parser.add_argument(
    '--position', dest='position', type=str,
    help='only examines the relevant terminal exon',
    default='all'
)
parser.add_argument(
    '--align', dest='align',
    help='Align highest points when calculating metaplot',
    action='store_true'
)

args = parser.parse_args()

def which(x,value=True):
    """Returns a list of locations in x that satisfy value"""
    return [a for a,b in enumerate(x) if b==value]


def notwhich(x,value=0):
    """Returns a list of locations in x that do not satisty value"""
    return [a for a,b in enumerate(x) if b!=value]


def flatten(list_of_lists):
    """Collapses a list/tuple of lists into a single list"""
    return [item for sublist in list_of_lists for item in sublist]


# Store bedgraph files in dicts
pos_coverage = {}
neg_coverage = {}

print('# Positive coverage data: {}'.format(args.pos_bg))
coverage_file = open(args.pos_bg)
for line in coverage_file:
    chrom,start,end,count = line.rstrip().split()
    if chrom not in pos_coverage:
        pos_coverage[chrom] = {}
    
    count = float(count)
    for i in range(int(start),int(end)):
        pos_coverage[chrom][i] = count

print('# Negative coverage data: {}'.format(args.neg_bg))
coverage_file = open(args.neg_bg)
for line in coverage_file:
    chrom,start,end,count = line.rstrip().split()
    if chrom not in neg_coverage:
        neg_coverage[chrom] = {}
    
    count = float(count)
    for i in range(int(start),int(end)):
        neg_coverage[chrom][i] = count

input_filetype = args.reference_GFF.split('.')[-1].lower()
if input_filetype == 'gff':
    input_filetype = 'gff3'

ref_transcripts = gu.get_file_content(args.reference_GFF,input_filetype)

ref_IDs = sorted(list(ref_transcripts.keys()))
print('# {} reference transcripts: {}'.format(
    len(ref_IDs),
    args.reference_GFF
))


ratio = 0
mratio = 0
ID_count = 0

metalength = args.flanking*2 + 1
start_meta_pos = [0]*metalength
start_meta_neg = [0]*metalength
end_meta_pos = [0]*metalength
end_meta_neg = [0]*metalength

for ID in ref_IDs:
    # current_transcript = ref_transcripts[ID].data
    current_transcript = ref_transcripts[ID]
    chrom = current_transcript.chrom
    strand = current_transcript.strand
    if strand != args.strand:
        continue
    
    if chrom not in pos_coverage:
        pos_coverage[chrom] = {}
        
    if chrom not in neg_coverage:
        neg_coverage[chrom] = {}
    
    # Gather information about exons
    # subfeatures = ref_transcripts[ID].subfeatures
    # exon_starts = sorted([int(i.data['start']) for i in subfeatures])
    # exon_ends = sorted([int(i.data['end']) for i in subfeatures])
    exon_starts = current_transcript.get_exon_start()
    exon_ends = current_transcript.get_exon_end()
    exon_starts[0] = exon_starts[0] - args.flanking
    exon_ends[-1] = exon_ends[-1] + args.flanking
    if exon_starts[0] < 0:
        exon_starts[0] = 0
    
    if args.position in ['5','5p','start','TSS']:
        # Check the strand of the transcript, and keep only the relevant terminal exon
        if strand == '+':
            positions = list(range(exon_starts[0] - 1, exon_ends[0]))
        elif strand == '-':
            positions = list(range(exon_starts[-1] - 1, exon_ends[-1]))
    elif args.position in ['3','3p','end','TES']:
        if strand == '+':
            positions = list(range(exon_starts[-1] - 1, exon_ends[-1]))
        elif strand == '-':
            positions = list(range(exon_starts[0] - 1, exon_ends[0]))
    else:
        # Make a list of all nucleotide positions in the input transcript
        positions = flatten(
            [list(range(a - 1, b)) for a,b in zip(exon_starts, exon_ends)]
        )
    
    positions = [i for i in positions if i > 0]
    
    # Get read coverage at each position with coverage
    pos = [
        pos_coverage[chrom].get(i,0)
        for i in positions
    ]
    
    neg = [
        neg_coverage[chrom].get(i,0)
        for i in positions
    ]
    
    if strand == '-':
        pos = list(reversed(pos))
        neg = list(reversed(neg))
    
    max_pos = max(pos)
    max_neg = max(neg)
    max_both = max([max_pos,max_neg])
    
    if max_pos < 1 or max_neg < 1:
        continue
    
    if args.align:
        # Offsets the center of the region to add to
        # the metaplot by its local maximum, resulting
        # in a sharper peak that reflects the signal spread
        # at each site rather than at sites generally
        start_window = pos[:metalength]
        if max(start_window) == 0:
            start_offset = 0
        else:
            start_offset = min(
                which(
                    start_window, 
                    max(start_window)
                )
            ) - args.flanking
        # Apply offset, adding zeros to pad
        # missing data when necessary
        if start_offset == 0:
            s_p_to_add = pos[:metalength]
            s_n_to_add = neg[:metalength]
        elif start_offset > 0:
            s_p_to_add = pos[start_offset:(metalength + start_offset)]
            s_n_to_add = neg[start_offset:(metalength + start_offset)]
            if len(s_p_to_add) < metalength:
                s_p_to_add += [0]*(metalength - len(s_p_to_add))
                s_n_to_add += [0]*(metalength - len(s_n_to_add))
        else:
            s_p_to_add = [0]*abs(start_offset) + pos[:(metalength + start_offset)]
            s_n_to_add = [0]*abs(start_offset) + neg[:(metalength + start_offset)]
            
        
        end_window = pos[-metalength:]
        if max(end_window) == 0:
            end_offset = 0
        else:
            end_offset = max(
                which(
                    end_window,
                    max(end_window)
                )
            ) - args.flanking
        # Apply offset, adding zeros to pad
        # missing data when necessary
        if end_offset == 0:
            e_p_to_add = pos[-metalength:]
            e_n_to_add = neg[-metalength:]
        elif end_offset > 0:
            e_p_to_add = pos[(end_offset - metalength):] + [0]*end_offset
            e_n_to_add = neg[(end_offset - metalength):] + [0]*end_offset
        else:
            e_p_to_add = pos[(end_offset - metalength):end_offset]
            e_n_to_add = neg[(end_offset - metalength):end_offset]
            if len(e_p_to_add) < metalength:
                e_p_to_add += [0]*(metalength - len(e_p_to_add))
                e_n_to_add += [0]*(metalength - len(e_n_to_add))
        
    else:
        s_p_to_add = pos[:metalength]
        s_n_to_add = neg[:metalength]
        e_p_to_add = pos[-metalength:]
        e_n_to_add = neg[-metalength:]

    # Add values from the current feature to the metaplots
    start_meta_pos = [
        a + (float(b)/max_both)
        for a,b in zip(start_meta_pos, s_p_to_add)
    ]
    start_meta_neg = [
        a + (float(b)/max_both)
        for a,b in zip(start_meta_neg, s_n_to_add)
    ]
    end_meta_pos = [
        a + (float(b)/max_both)
        for a,b in zip(end_meta_pos, e_p_to_add)
    ]
    end_meta_neg = [
        a + (float(b)/max_both)
        for a,b in zip(end_meta_neg, e_n_to_add)
    ]
    
    ratio += float(max_neg)/max_pos
    # mratio += float(mean_neg)/mean_pos
    ID_count += 1

print('# {} {} stranded transcripts with positive signal'.format(
    ID_count,
    args.strand
))

scaling_factor = ratio/ID_count

start_meta = [(a*scaling_factor - b)/ID_count for a,b in zip(start_meta_pos, start_meta_neg)]
end_meta = [(a*scaling_factor - b)/ID_count for a,b in zip(end_meta_pos, end_meta_neg)]

start_sum = sum([i for i in start_meta if i > 0])
if start_sum == 0:
    start_band = args.maxband
else:
    start_proportions = [sum([j for j in start_meta[(args.flanking-i):(args.flanking+i+1)] if j>0])/start_sum for i in range(args.maxband)]+[1]
    start_band = min(which([i>=.6827 for i in start_proportions]))

end_sum = sum([i for i in end_meta if i > 0])
if end_sum == 0:
    end_band = args.maxband
else:
    end_proportions = [sum([j for j in end_meta[(args.flanking-i):(args.flanking+i+1)] if j>0])/end_sum for i in range(args.maxband)]+[1]
    end_band = min(which([i>=.6827 for i in end_proportions]))

if start_band == 0:
    start_band = args.maxband

if end_band == 0:
    end_band = args.maxband
    
print('# scaling_factor\tstart_band\tend_band')
print('{}\t{}\t{}'.format(
    scaling_factor,
    start_band,
    end_band
    )
)
