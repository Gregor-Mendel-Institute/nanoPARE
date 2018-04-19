import os
import re
import sys
import argparse
import time
import fasta_utils as fu
import multiprocessing as mp

###################
# INPUT ARGUMENTS #
###################

parser = argparse.ArgumentParser()
parser.add_argument(
    "-S", "--samplename", dest='SAMPLENAME',
    help="Name of the sample, used as prefix for output files."
)
parser.add_argument(
    "-I", "--input", dest='SAMFILES', nargs='+',
    help="Path(s) to aligned SAM file(s). Require header."
)
parser.add_argument(
    "-R", "--readtype", dest='READTYPE',
    help="Type of reads in input SAM file(s)",
    default="BODY", type=str
)
parser.add_argument(
    "-E", "--endtype", dest='ENDTYPE',
    help="Which end of each read to store (5p or 3p)",
    default="5p", type=str, choices=['5p','3p']
)
parser.add_argument(
    "-F", "--fasta", dest='FASTA',
    help="Genome FASTA file",
    default=None, type=str
)
parser.add_argument(
    "--secondary", dest='SECONDARY',
    help="(True/False) Include secondary alignments.",
    default=False, action="store_true"
)
parser.add_argument(
    "--allow_nonstranded", dest='ALLOW_NS',
    help="Allow rescue based on nonstranded coverage.",
    default=False,action="store_true"
)
parser.add_argument(
    "--allow_naive", dest='ALLOW_NAIVE',
    help="Allow naive mapping (equal weight)",
    default=False, action="store_true"
)
parser.add_argument(
    "--parallel", dest='PARALLEL',
    help="Utilize multiprocessing to write bedgraph",
    default=False, action="store_true"
)
parser.add_argument(
    "--map_softclip", dest='MAP_SOFTCLIP',
    help="If true, maps 5' end regardless of softclipping",
    default=False, action="store_true"
)
parser.add_argument(
    "--bed", dest='WRITE_BED',
    help="If true, outputs a BED file for each unique read",
    default=False, action="store_true"
)
parser.add_argument(
    "--softclip_type", dest='SOFTCLIP_TYPE',
    help="Choose which softclipped reads to keep.",
    default='both', choices=['none','5p','3p','both']
)
parser.add_argument(
    "--mmIO", dest='mmIO', type=str,
    help="(out) write multimappers to disk. (in) import multimappers from current directory",
    default=False, choices=['in','out','False']
)
parser.add_argument(
    "--sj_out", dest='SJ_OUT', action='store_true',
    help="If true, outputs splice junction table",
    default=False
)
parser.add_argument(
    "--weighted", dest='WEIGHTED', action='store_true',
    help="If true, adjusts reads with weights from XW tag",
    default=False
)
parser.add_argument(
    "--size_classes", dest='SIZE_CLASSES',
    help="(optional) outputs a table of read lengths per genome position 5' end",
    default=[], nargs='+'
)
parser.add_argument(
    "--untemp_out", dest='UNTEMP_OUT',
    help="(optional) outputs a bedgraph of per-position untemplated nucleotides (choose ACGTN)",
    default=[], nargs='+'
)
parser.add_argument(
    "--minmatch", dest='MINMATCH',
    help="Minimum number of matching nucleotides for a read.",
    default=16, type=int
)

args = parser.parse_args()

######################################
# ENVIRONMENT SETUP: DATA STRUCTURES #
######################################

if args.SIZE_CLASSES:
    args.SIZE_CLASSES = sorted([int(i) for i in args.SIZE_CLASSES])
    
if args.UNTEMP_OUT:
    args.UNTEMP_OUT = sorted([i.upper() for i in args.UNTEMP_OUT])

if args.FASTA:
    print('Loading genome FASTA file...')
    genome = fu.import_genome(args.FASTA)

if args.UNTEMP_OUT and not args.FASTA:
    print('ERROR: Untemplated nucleotide analysis requires a reference genome. Use the --fasta argument.')
    sys.exit(1)

print('Loading chromosome lengths from header: ' + args.SAMFILES[0])
chromosomes = {}

file = open(args.SAMFILES[0])
line = file.readline()
while line[0] == '@':
    l = line.rstrip().split('\t')
    if l[0] == '@SQ':
        chrom = l[1][3:]
        chromlen = int(l[2][3:])
        chromosomes[chrom] = chromlen
    
    line = file.readline()
    if not line:
        break
    

file.close()
genome_size = sum(chromosomes.values())

# The dicts COVERAGE and ENDPOINT will store read count information
# to be written as BEDGRAPH files.
print('Constructing coverage dictionaries...')
COVERAGE = {}
ENDPOINT = {}
READLENGTHS = {}
UNTEMP = {}
SJ = {}

COVERAGE['+'] = {}
COVERAGE['-'] = {}

ENDPOINT['+'] = {}
ENDPOINT['-'] = {}

READLENGTHS['+'] = {}
READLENGTHS['-'] = {}

UNTEMP['+'] = {}
UNTEMP['-'] = {}

SJ['+'] = {}
SJ['-'] = {}
SJ['.'] = {}

for chrom,length in chromosomes.items():
    COVERAGE['+'][chrom] = {}
    COVERAGE['-'][chrom] = {}
    ENDPOINT['+'][chrom] = {}
    ENDPOINT['-'][chrom] = {}
    READLENGTHS['+'][chrom] = {}
    READLENGTHS['-'][chrom] = {}
    UNTEMP['+'][chrom] = {}
    UNTEMP['-'][chrom] = {}
    SJ['+'][chrom] = {}
    SJ['-'][chrom] = {}
    SJ['.'][chrom] = {}

# One class is defined for this script:
# RNASeqRead objects contain all information for a single read / read pair
# Mappings are stored as a quad in one of the lists 'primary' or 'secondary':
#   (pos, chrom, strand, cigar)

class RNASeqRead():
    def __init__(self,ID=None,paired=False,Nmap=0,primary={},secondary={},seq1=None,seq2=None):
        self.paired = paired
        self.Nmap = Nmap
        self.primary = primary
        self.secondary = secondary
        self.ID = ID
        self.seq1 = seq1
        self.seq2 = seq2
    
    def __repr__(self):
        return '<RNASeqRead object>'
    
    def __eq__(self, other) : 
        return self.__dict__ == other.__dict__
    
    def map_locations(self,include_secondary=args.SECONDARY):
        mappings = list(self.primary.values())
        if include_secondary:
            if self.secondary:
                mappings += list(self.secondary.values())
        return mappings
        
    def minpos(self,include_secondary=args.SECONDARY):
        positions = list(self.primary.keys())
        if include_secondary:
            if self.secondary:
                positions.extend(list(self.secondary.keys()))
        return min(positions)
        
    def mincigar(self,include_secondary=args.SECONDARY):
        minpos = self.minpos(include_secondary)
        map = self.primary.get(minpos,None)
        if not map:
            map = self.secondary.get(minpos,None)
        return sorted(map[2])[0][1]

# The dict MM retains all MultiMappers from the input libraries,
# subclassifying them based on the number of times they mapped.
# This allows for rapid hierarchical recall of multimapping reads
# from fewest to most mapping locations.

MM = {}

################################
# ENVIRONMENT SETUP: FUNCTIONS #
################################
    
def flatten(list_of_lists):
    return [item for sublist in list_of_lists for item in sublist]

def which(x,value=True):
    return [a for a,b in enumerate(x) if b==value]

    
def writer(merge_filename,shared_queue,stop_token):
    """Initializes a writer to output the mp queue."""
    dest_file = open(merge_filename,'w')
    while True:
        line = shared_queue.get()
        if line == stop_token:
            dest_file.close()
            return
        dest_file.write(line)


def get_junction_strand(junction,chrom):
    """ Returns '+', '-', or '.' for a left/right
    pair of splice junction positions based
    on the flanking genomic sequence """
    left,right = junction
    flanking_sequence = genome[chrom][left-1:left+1] + genome[chrom][right-2:right]
    flanking_sequence = flanking_sequence.upper()
    if flanking_sequence in ['GTAG','GCAG','ATAC']:
        strand = '+'
    elif flanking_sequence in ['CTAC','CTGC','GTAT']:
        strand = '-'
    elif flanking_sequence == 'GTAC':
        strand = '.'
    # elif flanking_sequence[:2] == 'GT' or flanking_sequence[-2:] == 'AG':
        # strand = '+'
    # elif flanking_sequence[:2] == 'CT' or flanking_sequence[-2:] == 'AC':
        # strand = '-'
    else:
        strand = '.'
    return strand
    
def populate_bedgraphs(read_object,mmIO=args.mmIO):
    """Adds values from read_object to the appropriate bedgraph dicts."""
    if len(read_object.primary)==0:
        return
    
    # If multimapper, store in MM and do not change bedgraph dicts
    nmap = len(read_object.map_locations())
    if nmap > 1:
        if mmIO == "out":
            return nmap
        
        minpos = read_object.minpos()
        cigar = read_object.mincigar()
        seq = None
        if read_object.seq1:
            seq = read_object.seq1
        if read_object.seq2:
            seq += read_object.seq2
        if seq:
            identifier = seq
        else:
            identifier = cigar
        
        read_object.ID = None
        if nmap in MM:
            if minpos in MM[nmap]:
                if identifier in MM[nmap][minpos]:
                    MM[nmap][minpos][identifier] = (read_object,MM[nmap][minpos][identifier][1]+1)
                else:
                    MM[nmap][minpos][identifier] = (read_object,1)
            else:
                MM[nmap][minpos] = {identifier:(read_object,1)}
        else:
            MM[nmap] = {minpos:{identifier:(read_object,1)}}
    else:
        # Parse the mapping information of the uniquely mapping read
        chrom,strand,poslist = read_object.map_locations()[0]
        map_positions = set()
        readweight = float(0)
        for pos,cigar,mate,weight in poslist:
            mate_positions,junctions = get_mapping_positions(pos,cigar)
            readweight += weight
            if junctions:
                # TODO: Add to a splice junction dictionary
                for j in junctions:
                    jstrand = get_junction_strand(j,chrom)
                    SJ[jstrand][chrom][j] = SJ[jstrand][chrom].get(j,0) + 1
            
            untemp_nucs = {}
            if args.UNTEMP_OUT:
                # If untemplated nucleotides are examined...
                # 1) Collect all softclipped nucleotides
                # 2) Does the softclipped nucleotide match the genome?
                #    Y: Is it contiguous with mate_positions?
                #       Y: Add to mate_positions
                #       N: Keep as untemplated
                #    N: Keep as untemplated
                untemp_positions = set()
                if strand == '+':
                    if mate == 1:
                        softclip_positions = get_untemp_positions(
                            pos,
                            cigar,
                            seq=read_object.seq1
                        )
                    elif mate == 2:
                        softclip_positions = get_untemp_positions(
                            pos,
                            cigar,
                            seq=fu.rc(read_object.seq2)
                        )
                        softclip_positions = dict(
                            [
                            (k,fu.IUPACcomp[v])
                            for k,v in softclip_positions.items()
                            ]
                        )
                    else:
                        print("ERROR: mate pair incorrectly specified.")
                        sys.exit(1)
                elif strand == '-':
                    if mate == 1:
                        softclip_positions = get_untemp_positions(
                            pos,
                            cigar,
                            seq=fu.rc(read_object.seq1)
                        )
                        softclip_positions = dict(
                            [
                            (k,fu.IUPACcomp[v])
                            for k,v in softclip_positions.items()
                            ]
                        )
                    elif mate == 2:
                        softclip_positions = get_untemp_positions(
                            pos,
                            cigar,
                            seq=read_object.seq2
                        )
                    else:
                        print("ERROR: mate pair incorrectly specified.")
                        sys.exit(1)
                    
                
                if softclip_positions:
                    if all([i > max(mate_positions) for i in softclip_positions.keys()]):
                        # Tail clip, sort from bottom to top
                        sp = sorted([(k,v) for k,v in softclip_positions.items()],reverse=False)
                    else:
                        # Head clip or mixed clip, sort from top to bottom
                        sp = sorted([(k,v) for k,v in softclip_positions.items()],reverse=True)
                    
                    for p,n in sp:
                        if p-1 < 0 or p-1 >= chromosomes.get(chrom,0):
                            continue
                        if chrom not in genome:
                            print('WARNING: chromosome {} not found'.format(chrom))
                            continue
                        if strand == '-':
                            nuc = fu.IUPACcomp[n]
                        else:
                            nuc = n
                        
                        if genome[chrom][p-1] == nuc:
                            # Softclipped nucleotide matches the genome
                            if p - 1 in mate_positions or p + 1 in mate_positions:
                                mate_positions.add(p)
                                # print('Templated {} found: {}:{} {}'.format(n,chrom,p,strand))
                            else:
                                untemp_positions.add(p)
                        else:
                            untemp_positions.add(p)
                    
                    if untemp_positions:
                        untemp_nucs = dict(
                            [
                                (k,v)
                                for k,v in softclip_positions.items()
                                if k in untemp_positions
                            ]
                        )
                        
                        if args.SOFTCLIP_TYPE == 'none':
                            # If any softclipping occurred, throw out the whole read
                            mate_positions = set()
                            continue
                        
                        sc_right = all([i > max(mate_positions) for i in untemp_nucs.keys()])
                        sc_left = all([i < min(mate_positions) for i in untemp_nucs.keys()])
                        if args.SOFTCLIP_TYPE == '5p':
                            # Look for 3p softclipping and discard
                            if strand == '+' and sc_right:
                                mate_positions = set()
                                continue
                            elif strand == '-' and sc_left:
                                mate_positions = set()
                                continue
                        elif args.SOFTCLIP_TYPE == '3p':
                            # Look for 5p softclipping and discard
                            if strand == '+' and sc_left:
                                mate_positions = set()
                                continue
                            elif strand == '-' and sc_right:
                                mate_positions = set()
                                continue
                        
                        if not all([v in args.UNTEMP_OUT for k,v in untemp_nucs.items()]):
                            # If the untemplated nucleotides aren't those being examined, discard the read
                            mate_positions = set()
                            continue
                        
                        # Update the UNTEMP dict with the correct location and identity of untemplated nucleotides
                        for p,n in untemp_nucs.items():
                            u_pos = p - 1
                            if u_pos not in UNTEMP[strand][chrom]:
                                UNTEMP[strand][chrom][u_pos] = {}
                            UNTEMP[strand][chrom][u_pos][n] = UNTEMP[strand][chrom][u_pos].get(n,0) + 1
                    
            mate_positions = set([i for i in mate_positions if i > 0 and i <= chromosomes[chrom]])
            map_positions.update(mate_positions)
        
        readweight = readweight/len(poslist)
        if len(map_positions) > args.MINMATCH:
            if strand == '+':
                if args.ENDTYPE == '5p':
                    ENDPOINT_end = min(map_positions) - 1
                else:
                    ENDPOINT_end = max(map_positions) - 1
            elif strand == '-':
                if args.ENDTYPE == '5p':
                    ENDPOINT_end = max(map_positions) - 1
                else:
                    ENDPOINT_end = min(map_positions) - 1
            
            
            if args.WRITE_BED:
                # Generate a line for a temp BED file
                temp_seq = ''
                temp_positions = sorted(list(map_positions))
                for i in range(len(temp_positions)):
                    if strand == '+':
                        nuc_to_add = genome[chrom][temp_positions[i] - 1]
                    elif strand == '-':
                        nuc_to_add = fu.IUPACcomp[genome[chrom][temp_positions[i] - 1]]
                    if i > 0:
                        if temp_positions[i] > 1 + temp_positions[i-1]:
                            # non-adjacent nucleotides are separated by a dot
                            temp_seq += '.'
                    temp_seq += nuc_to_add
                
                left_seq = ''
                right_seq = ''
                if args.UNTEMP_OUT and untemp_nucs:
                    left_nucs = dict([(p,n) for p,n in untemp_nucs.items() if p < min(temp_positions)])
                    right_nucs = dict([(p,n) for p,n in untemp_nucs.items() if p > max(temp_positions)])
                    if left_nucs:
                        left_positions = sorted(list(left_nucs.keys()))
                        for i in range(len(left_positions)):
                            nuc_to_add = left_nucs[left_positions[i]]
                            if i > 0:
                                if left_positions[i] > 1 + left_positions[i-1]:
                                    # non-adjacent nucleotides are separated by a dot
                                    left_seq += '.'
                            left_seq += nuc_to_add
                    if right_nucs:
                        right_positions = sorted(list(right_nucs.keys()))
                        for i in range(len(right_positions)):
                            nuc_to_add = right_nucs[right_positions[i]]
                            if i > 0:
                                if right_positions[i] > 1 + right_positions[i-1]:
                                    # non-adjacent nucleotides are separated by a dot
                                    right_seq += '.'
                            right_seq += nuc_to_add
                        
                if strand == '+':
                    head_seq = left_seq
                    tail_seq = right_seq
                elif strand == '-':
                    # Reverse the templated string
                    head_seq = right_seq[::-1]
                    temp_seq = temp_seq[::-1]
                    tail_seq = left_seq[::-1]
                
                full_read = head_seq + temp_seq + tail_seq
                read_length = len(''.join(full_read.split('.')))
                
                # chrom	start	end	unique_name	read_num	strand	head [untemp 5']	templated sequence	tail  [untemp 3']
                temp_bed.write(
                    '\t'.join(
                        [
                            str(i) for i in [
                                chrom,
                                min(temp_positions) - 1,
                                max(temp_positions),
                                'tmp',
                                float(1)*readweight,
                                strand,
                                head_seq,
                                temp_seq,
                                tail_seq,
                                len(read_object.seq1)
                            ]
                        ]
                    ) + '\n'
                )
            
            ENDPOINT[strand][chrom][ENDPOINT_end] = ENDPOINT[strand][chrom].get(ENDPOINT_end,0) + float(1)*readweight
            if args.SIZE_CLASSES:
                n = len(read_object.seq1)
                if ENDPOINT_end not in READLENGTHS[strand][chrom]:
                    READLENGTHS[strand][chrom][ENDPOINT_end] = {}
                READLENGTHS[strand][chrom][ENDPOINT_end][n] = READLENGTHS[strand][chrom][ENDPOINT_end].get(n,0) + float(1)*readweight
            for pos in map_positions:
                COVERAGE[strand][chrom][pos-1] = COVERAGE[strand][chrom].get(pos-1,0) + float(1)*readweight
        if mmIO == "out":
            return 0


def assign_multimapper(
    read_object,
    value=1,
    allow_naive=False,
    include_secondary=args.SECONDARY,
    allow_nonstranded=args.ALLOW_NS
):
    """Assigns values to the respective bedgraphs for reads mapping to multiple locations
    using a 'rich-get-richer' algorithm:
    1) If stranded data of the same read type exists within region of coverage, assign proportionally
    2) Else if coverage data of any type exists, assign proportionally
    3) Else if allow_naive, then distribute evenly across all mapping positions
    4) Else do nothing, returning the original read_object
    """
    # Get dictionary of all mappings from the read_object
    all_mappings = read_object.map_locations()
    
    mappos_list = []
    junction_list = []
    weights_list = []
    for chrom,strand,poslist in all_mappings:
        allpositions = set()
        alljunctions = []
        currentweight = 0
        for p,c,m,w in poslist:
            positions, junctions = get_mapping_positions(p,c)
            allpositions.update(positions)
            currentweight += w
            for j in junctions:
                if j not in alljunctions:
                    alljunctions += [j]
        
        weights_list += [currentweight/len(poslist)]
        mappos_list += [allpositions]
        junction_list += [alljunctions]
    
    
    long_enough = [len(i) >= args.MINMATCH for i in mappos_list]
    all_mappings = [all_mappings[i] for i in which(long_enough)]
    
    mappos_list = [mappos_list[i] for i in which(long_enough)]
    junction_list = [junction_list[i] for i in which(long_enough)]
    weights_list = [weights_list[i] for i in which(long_enough)]
    chroms_list = [chrom for chrom,strand,poslist in all_mappings]
    strand_list = [strand for chrom,strand,poslist in all_mappings]
    poslist_list = [poslist for chrom,strand,poslist in all_mappings]
    
    if len(all_mappings) == 0:
        return
    
    # Determine existing reads
    existing_reads = [
        sum(
            [
                COVERAGE[strand_list[i]][chroms_list[i]].get(p-1,0)
                for p in mappos_list[i]
            ]
        ) for i in range(len(all_mappings))
    ]
    
    if allow_nonstranded:
        if sum(existing_reads) == 0:
            oppstrands = {'+':'-','-':'+'}
            oppstrand_list = [oppstrands[i] for i in strand_list]
            existing_reads = [
                sum(
                    [
                        COVERAGE[oppstrand_list[i]][chroms_list[i]].get(p-1,0)
                        for p in mappos_list[i]
                    ]
                ) for i in range(len(all_mappings))
            ]
    
    if sum(existing_reads) == 0:
        if allow_naive:
            existing_reads = [1]*len(all_mappings)
        else:
            return (read_object,value)

    # Assign fractional reads
    total_coverage = float(sum(existing_reads))
    proportions = [float(i)/total_coverage for i in existing_reads]
    untemp_dict = {}
    
    for i in range(len(all_mappings)):

        current_strand = strand_list[i]
        
        if args.UNTEMP_OUT:
            # If untemplated nucleotides are examined...
            # 1) Collect all softclipped nucleotides
            # 2) Does the softclipped nucleotide match the genome?
            #    Y: Is it contiguous with mappos_list[i]?
            #       Y: Add to mappos_list[i]
            #       N: Keep as untemplated
            #    N: Keep as untemplated
            for m in range(len(poslist_list[i])):
                untemp_positions = set()
                pos,cigar,mate,weight = poslist_list[i][m]
                softclip_positions = None
                if current_strand == '+':
                    if mate == 1:
                        softclip_positions = get_untemp_positions(
                            pos,
                            cigar,
                            seq=read_object.seq1
                        )
                    elif mate == 2:
                        softclip_positions = get_untemp_positions(
                            pos,
                            cigar,
                            seq=fu.rc(read_object.seq2)
                        )
                        softclip_positions = dict(
                            [
                            (k,fu.IUPACcomp[v])
                            for k,v in softclip_positions.items()
                            ]
                        )
                    else:
                        print("ERROR: mate pair incorrectly specified.")
                        sys.exit(1)
                elif current_strand == '-':
                    if mate == 1:
                        softclip_positions = get_untemp_positions(
                            pos,
                            cigar,
                            seq=fu.rc(read_object.seq1)
                        )
                        softclip_positions = dict(
                            [
                            (k,fu.IUPACcomp[v])
                            for k,v in softclip_positions.items()
                            ]
                        )
                    elif mate == 2:
                        softclip_positions = get_untemp_positions(
                            pos,
                            cigar,
                            seq=read_object.seq2
                        )
                    else:
                        print("ERROR: mate pair incorrectly specified.")
                        sys.exit(1)
                    
                
                if softclip_positions:
                    if all([k > max(mappos_list[i]) for k in softclip_positions.keys()]):
                        # Tail clip, sort from bottom to top
                        sp = sorted([(k,v) for k,v in softclip_positions.items()],reverse=False)
                    else:
                        # Head clip or mixed clip, sort from top to bottom
                        sp = sorted([(k,v) for k,v in softclip_positions.items()],reverse=True)
                    
                    for p,n in sp:
                        if p-1 < 0 or p-1 >= chromosomes.get(chroms_list[i],0):
                            continue
                        if chroms_list[i] not in genome:
                            print('WARNING: chromosome {} not found'.format(chroms_list[i]))
                            continue
                        if current_strand == '-':
                            nuc = fu.IUPACcomp[n]
                        else:
                            nuc = n
                        
                        if genome[chroms_list[i]][p-1] == nuc:
                            # Softclipped nucleotide matches the genome
                            if p - 1 in mappos_list[i] or p + 1 in mappos_list[i]:
                                mappos_list[i].add(p)
                            else:
                                untemp_positions.add(p)
                        else:
                            untemp_positions.add(p)
                    
                    if untemp_positions:
                        untemp_nucs = dict(
                            [
                                (k,v)
                                for k,v in softclip_positions.items()
                                if k in untemp_positions
                            ]
                        )
                        
                        if args.SOFTCLIP_TYPE == 'none':
                            # If any softclipping occurred, throw out the whole read
                            proportions[i] = 0
                            continue
                        
                        sc_right = all([k > max(mappos_list[i]) for k in untemp_nucs.keys()])
                        sc_left = all([k < min(mappos_list[i]) for k in untemp_nucs.keys()])
                        if args.SOFTCLIP_TYPE == '5p':
                            # Look for 3p softclipping and discard
                            if current_strand == '+' and sc_right:
                                proportions[i] = 0
                                continue
                            elif current_strand == '-' and sc_left:
                                proportions[i] = 0
                                continue
                        elif args.SOFTCLIP_TYPE == '3p':
                            # Look for 5p softclipping and discard
                            if current_strand == '+' and sc_left:
                                proportions[i] = 0
                                continue
                            elif current_strand == '-' and sc_right:
                                proportions[i] = 0
                                continue
                        
                        if not all([v in args.UNTEMP_OUT for k,v in untemp_nucs.items()]):
                            # If the untemplated nucleotides aren't those being examined, discard the read
                            proportions[i] = 0
                            continue
                        
                        untemp_dict[i] = untemp_nucs
    
    # Readjust proportions in case some mappings were removed
    if sum(proportions) == 0:
        return
    proportions = [float(j)/sum(proportions) for j in proportions]
    
    for i in range(len(all_mappings)):
        if proportions[i] == 0:
            continue
        
        if strand_list[i] == '+':
            if args.ENDTYPE == '5p':
                ENDPOINT_end = min(mappos_list[i]) - 1
            else:
                ENDPOINT_end = max(mappos_list[i]) - 1
        elif strand_list[i] == '-':
            if args.ENDTYPE == '5p':
                ENDPOINT_end = max(mappos_list[i]) - 1
            else:
                ENDPOINT_end = min(mappos_list[i]) - 1
        
        if args.WRITE_BED:
            # Write a line to a temp BED file
            temp_seq = ''
            temp_positions = sorted(list(mappos_list[i]))
            for p in range(len(temp_positions)):
                if strand_list[i] == '+':
                    nuc_to_add = genome[chroms_list[i]][temp_positions[p] - 1]
                elif strand_list[i] == '-':
                    nuc_to_add = fu.IUPACcomp[genome[chroms_list[i]][temp_positions[p] - 1]]
                if p > 0:
                    if temp_positions[p] > 1 + temp_positions[p-1]:
                        # non-adjacent nucleotides are separated by a dot
                        temp_seq += '.'
                temp_seq += nuc_to_add
            
            left_seq = ''
            right_seq = ''
            if args.UNTEMP_OUT and i in untemp_dict:
                left_nucs = dict([(p,n) for p,n in untemp_dict[i].items() if p < min(temp_positions)])
                right_nucs = dict([(p,n) for p,n in untemp_dict[i].items() if p > max(temp_positions)])
                if left_nucs:
                    left_positions = sorted(list(left_nucs.keys()))
                    for p in range(len(left_positions)):
                        nuc_to_add = left_nucs[left_positions[p]]
                        if p > 0:
                            if left_positions[p] > 1 + left_positions[p-1]:
                                # non-adjacent nucleotides are separated by a dot
                                left_seq += '.'
                        left_seq += nuc_to_add
                if right_nucs:
                    right_positions = sorted(list(right_nucs.keys()))
                    for p in range(len(right_positions)):
                        nuc_to_add = right_nucs[right_positions[p]]
                        if p > 0:
                            if right_positions[p] > 1 + right_positions[p-1]:
                                # non-adjacent nucleotides are separated by a dot
                                right_seq += '.'
                        right_seq += nuc_to_add
                    
            if strand_list[i] == '+':
                head_seq = left_seq
                tail_seq = right_seq
            elif strand_list[i] == '-':
                # Reverse the templated string
                head_seq = right_seq[::-1]
                temp_seq = temp_seq[::-1]
                tail_seq = left_seq[::-1]
            
            full_read = head_seq + temp_seq + tail_seq
            read_length = len(''.join(full_read.split('.')))
            
            # chrom	start	end	unique_name	read_num	strand	head [untemp 5']	templated sequence	tail  [untemp 3']
            temp_bed.write(
                '\t'.join(
                    [
                        str(j) for j in [
                            chroms_list[i],
                            min(temp_positions) - 1,
                            max(temp_positions),
                            'tmp',
                            float(value)*proportions[i]*weights_list[i],
                            strand_list[i],
                            head_seq,
                            temp_seq,
                            tail_seq,
                            len(read_object.seq1)
                        ]
                    ]
                ) + '\n'
            )
        
        junctions = junction_list[i]
        if junctions:
            # TODO: Add to a splice junction dictionary
            for j in junctions:
                jstrand = get_junction_strand(j,chroms_list[i])
                SJ[jstrand][chrom][j] = SJ[jstrand][chrom].get(j,0) + \
                    float(value)*proportions[i]*weights_list[i]
        
        ENDPOINT[strand_list[i]][chroms_list[i]][ENDPOINT_end] = \
            ENDPOINT[strand_list[i]][chroms_list[i]].get(ENDPOINT_end,0) + \
            float(value)*proportions[i]*weights_list[i]
        
        if args.SIZE_CLASSES:
            n = len(read_object.seq1)
            if ENDPOINT_end not in READLENGTHS[strand_list[i]][chroms_list[i]]:
                READLENGTHS[strand_list[i]][chroms_list[i]][ENDPOINT_end] = {}
            READLENGTHS[strand_list[i]][chroms_list[i]][ENDPOINT_end][n] = \
                READLENGTHS[strand_list[i]][chroms_list[i]][ENDPOINT_end].get(n,0) + \
                float(value)*proportions[i]*weights_list[i]
        
        for pos in mappos_list[i]:
            COVERAGE[strand_list[i]][chroms_list[i]][pos-1] = \
                COVERAGE[strand_list[i]][chroms_list[i]].get(pos-1,0) + \
                float(value)*proportions[i]*weights_list[i]
        
        # Update the UNTEMP dict with the correct location and identity of untemplated nucleotides
        if i in untemp_dict:
            for p,n in untemp_dict[i].items():
                u_pos = p - 1
                if u_pos not in UNTEMP[strand_list[i]][chroms_list[i]]:
                    UNTEMP[strand_list[i]][chroms_list[i]][u_pos] = {}
                
                UNTEMP[strand_list[i]][chroms_list[i]][u_pos][n] = \
                    UNTEMP[strand_list[i]][chroms_list[i]][u_pos].get(n,0) + \
                    float(value)*proportions[i]*weights_list[i]


def run_mm_assignment(multiplicity):
    oscillator = False
    reads_were_assigned = True
    while reads_were_assigned:
        reads_were_assigned = False
        mm_readsort = sorted(
            list(MM[multiplicity].keys()),
            reverse=oscillator
        )
        keep = []
        for pos in mm_readsort:
            current_read = MM[multiplicity].get(pos,{})
            if not current_read:
                continue
            failed_to_map = {}
            for k in current_read.keys():
                read_object,value = current_read[k]
                unmapped = assign_multimapper(read_object,value,allow_naive=False)
                if unmapped:
                    failed_to_map[k] = unmapped
            MM[multiplicity][pos] = failed_to_map
            if failed_to_map:
                keep.append(pos)
            else:
                del MM[multiplicity][pos]
                reads_were_assigned = True
        oscillator = not oscillator
    
    # Assign any remaining reads in the current multiplicity group naively
    residuals = sum(
        [
            count for obj,count in flatten(
                [a.values() for a in MM[multiplicity].values()]
            )
        ]
    )
    
    if args.ALLOW_NAIVE:
        mm_readsort = sorted(list(MM[multiplicity].keys()))
        keep = []
        for pos in mm_readsort:
            current_read = MM[multiplicity].get(pos,{})
            if not current_read:
                continue
            failed_to_map = {}
            for k in current_read.keys():
                read_object,value = current_read[k]
                unmapped = assign_multimapper(read_object,value,allow_naive=True)
                if unmapped:
                    failed_to_map[k] = unmapped
            MM[multiplicity][pos] = failed_to_map
            if failed_to_map:
                keep.append(pos)
            else:
                del MM[multiplicity][pos]
                reads_were_assigned = True
        
        assert len(MM[multiplicity]) == 0, "ERROR: reads still remaining in {}:{}".format(args.READTYPE,multiplicity)
        
    return residuals


def write_bedgraph_from_dict(input, output_filename, digits=2, parallel=False, multi_key=False):
    """Writes unsorted BEDGRAPH to output_filename from input dict"""
    
    def writer(merge_filename,shared_queue,stop_token):
        """Initializes a writer to output the mp queue."""
        dest_file = open(merge_filename,'w')
        while True:
            line = shared_queue.get()
            if line == stop_token:
                dest_file.close()
                return
            dest_file.write(line)
    
    
    def generate_bedgraph_lines(values_dict, chrom, queue, digits=digits, parallel=parallel, multi_key=multi_key):
        """Converts dict of values to bedgraph format"""
        start = 0
        prevpos = 0
        prevcount = None
        chromlen = len(input[chrom])
        position = None
        # Iterate through every position with values in the dict
        for position in sorted(list(values_dict.keys())):
            if multi_key:
                all_values = [str(round(values_dict[position].get(k,0),digits)) for k in multi_key]
                count = '\t'.join(all_values)
            else:
                count = round(values_dict[position],digits)
            
            if count != prevcount or int(position) > 1 + prevpos:
                # The newly encountered value is not a continuation
                # of the previous value. Write the old run and start another.
                if prevcount and prevcount != 0:
                    line_to_write = '\t'.join(
                        [
                            str(i) for i in [
                                chrom,
                                start,
                                prevpos + 1,
                                prevcount
                            ]
                        ]
                    ) + '\n'
                    
                    if parallel:
                        # Write the old run to outfile
                        queue.put(line_to_write)
                    else:
                        # Write the old run to outfile
                        queue.write(line_to_write)
                        
                start = position
            prevcount = count
            prevpos = int(position)
        
        if position and prevcount and prevcount != 0:
            line_to_write = '\t'.join(
                [
                    str(i) for i in [
                        chrom,
                        start,
                        prevpos + 1,
                        prevcount
                    ]
                ]
            ) + '\n'
            
            if parallel:
                queue.put(line_to_write)
            else:
                queue.write(line_to_write)
    
    
    if parallel:
        queue = mp.Queue()
        STOP_TOKEN = "FINISHED"
        writer_process = mp.Process(
            target=writer,
            args=(output_filename,queue,STOP_TOKEN)
        )
        writer_process.start()
        
        if multi_key:
            first_line = '#chrom\tstart\tend\t{}\n'.format(
                '\t'.join([str(i) for i in multi_key])
            )
            queue.put(first_line)
        
        all_threads = []
        
        for chrom in sorted(list(input.keys())):
            all_threads.append(
                mp.Process(
                    target=generate_bedgraph_lines,
                    args=(
                        input[chrom],
                        chrom,
                        queue
                    )
                )
            )
        for i in range(len(all_threads)):
            all_threads[i].start()
        while len(mp.active_children()) > 1:
            time.sleep(1)
        queue.put("FINISHED")
        while len(mp.active_children()) > 0:
            time.sleep(1)
    else:
        queue = open(output_filename, 'w')
        if multi_key:
            first_line = '#chrom\tstart\tend\t{}\n'.format(
                '\t'.join([str(i) for i in multi_key])
            )
            queue.write(first_line)
        
        for chrom in sorted(list(input.keys())):
            generate_bedgraph_lines(
                input[chrom],
                chrom,
                queue,
                parallel = False
            )
        queue.close()


def get_mapping_positions(pos,cigar):
    """Converts the pos+CIGAR string to a set of mapped locations
    
    (description from http://bioinformatics.cvr.ac.uk/blog/tag/cigar-string/)
    CIGAR operators:
    D    Deletion; present in the reference but not in the read
    H    Hard Clipping; not present in the read
    I    Insertion; present in the read  but not in the reference
    M    Match; can be either an alignment match or mismatch
    N    Skipped region; a region is not present in the read
    P    Padding; padded area in the read and not in the reference
    S    Soft Clipping;  the clipped nucleotides are present in the read
    """
    positions = set()
    junctions = []
    current_pos = pos
    positions.add(current_pos)
    cigar_elements = zip(re.findall('[A-Z]+',cigar),[int(i) for i in re.findall('\d+',cigar)])
    first_element = True
    for operator,number in cigar_elements:
        if operator == 'S':
            if args.MAP_SOFTCLIP:
                if first_element:
                    positions.update(range(current_pos-number,current_pos))
                else:
                    positions.update(range(current_pos,current_pos+number))
            else:
                continue
        
        first_element = False
        if operator == 'M':
            positions.update(range(current_pos,current_pos+number))
            current_pos += number
        elif operator == 'N':
            if args.SJ_OUT:
                # TODO: Add donor and acceptor sites to read_object mapping
                leftside = current_pos
                current_pos += number
                rightside = current_pos - 1
                junctions += [(leftside,rightside)]
            else:
                current_pos += number
        elif operator == 'D':
            current_pos += number
        
        elif operator in ['I','H']:
            continue
    
    return positions, junctions
    
def get_untemp_positions(pos,cigar,seq,include_mismatch=False):
    """Converts the pos+CIGAR string to a dict of softclipped locations"""
    untemp_positions = {}
    # Track two positions:
    # genome_pos: coordinates of that nucleotide on a chromosome
    # string_pos: coordinates of the current nucleotide in the sequenced string
    genome_pos = pos
    string_pos = 0
    
    cigar_elements = zip(re.findall('[A-Z]+',cigar),[int(i) for i in re.findall('\d+',cigar)])
    first_element = True
    for operator,number in cigar_elements:
        if operator == 'S':
            if first_element:
                softclipped_positions = range(genome_pos-number,genome_pos)
                for p in range(len(softclipped_positions)):
                    untemp_positions[softclipped_positions[p]] = seq[string_pos + p]
            else:
                softclipped_positions = range(genome_pos,genome_pos+number)
                for p in range(len(softclipped_positions)):
                    untemp_positions[softclipped_positions[p]] = seq[string_pos + p]
            string_pos += number
        
        first_element = False
        if operator == 'M':
            genome_pos += number
            string_pos += number
        elif operator in ['D','N']:
            genome_pos += number
        elif operator in ['I']:
            string_pos += number
        elif operator in ['H']:
            continue
    
    return untemp_positions
        

def generate_read_from_sam(input_lines,readtype=args.READTYPE,keep_seq=False):
    """Convert a set of lines from a SAM file into an RNASeqRead object"""
    generated_read = RNASeqRead(ID=None,Nmap=0,primary={},secondary={},seq1=None,seq2=None)
    for line in input_lines:
        l = line.rstrip().split('\t')
        ID = l[0]
        if not generated_read.ID:
            generated_read.ID = ID
        assert ID == generated_read.ID, 'ERROR: nonmatching IDs in input_lines:\n{}\t{}'.format(generated_read.ID,ID)
        
        attributes = dict([(':'.join(i.split(':')[0:2]),i.split(':')[-1]) for i in l[11:]])
        if args.WEIGHTED:
            weight = float(attributes.get('XW:f',1))
        else:
            weight = float(1)
        
        # Nmap = int(attributes['NH:i'])
        # if generated_read.Nmap == 0:
            # generated_read.Nmap = Nmap
        # assert Nmap == generated_read.Nmap, 'ERROR: inconsistent Nmap score for {}'.format(ID)
        
        SAMflags   = bin(int(l[1]))[2:]
        SAMflags   = '0'*(12-len(SAMflags))+SAMflags
        # Interpret the binary SAM flags
        is_paired      = bool(int(SAMflags[-1]))
        pair_is_mapped = bool(int(SAMflags[-2]))
        read_reverse   = bool(int(SAMflags[-5]))
        first_in_pair  = bool(int(SAMflags[-7]))
        secondary      = bool(int(SAMflags[-9]))
        supplementary  = bool(int(SAMflags[-12]))
        # read_unmapped  = bool(int(SAMflags[-3]))
        # mate_unmapped  = bool(int(SAMflags[-4]))
        # mate_reverse   = bool(int(SAMflags[-6]))
        # second_in_pair = bool(int(SAMflags[-8]))
        # quality_fail   = bool(int(SAMflags[-10]))
        # pcr_duplicate  = bool(int(SAMflags[-11]))
        ###
        chrom      = l[2]
        pos        = int(l[3])
        pairpos    = int(l[7])
        mapscore   = int(l[4])
        cigar      = l[5]
        seq        = l[9]
        
        generated_read.paired = is_paired
        
        if keep_seq:
            # store the + stranded read of the respective mate pair
            if read_reverse:
                rcseq = fu.rc(seq)
                if first_in_pair or not is_paired:
                    # store reverse complement in seq1
                    if not generated_read.seq1:
                        generated_read.seq1 = rcseq
                    assert rcseq == generated_read.seq1, 'ERROR: nonmatching sequence in input_lines:\n{}\t{}'.format(generated_read.seq1,rcseq)
                else:
                    # store reverse complement in seq2
                    if not generated_read.seq2:
                        generated_read.seq2 = rcseq
                    assert rcseq == generated_read.seq2, 'ERROR: nonmatching sequence in input_lines:\n{}\t{}'.format(generated_read.seq2,rcseq)
            else:
                if first_in_pair or not is_paired:
                    # store in seq1
                    if not generated_read.seq1:
                        generated_read.seq1 = seq
                    assert seq == generated_read.seq1, 'ERROR: nonmatching sequence in input_lines:\n{}\t{}'.format(generated_read.seq1,seq)
                else:
                    # store in seq2
                    if not generated_read.seq2:
                        generated_read.seq2 = seq
                    assert seq == generated_read.seq2, 'ERROR: nonmatching sequence in input_lines:\n{}\t{}'.format(generated_read.seq2,seq)
        
        if read_reverse:
            if first_in_pair or not is_paired:
                strand = '-'
            else:
                strand = '+'
        else:
            if first_in_pair or not is_paired:
                strand = '+'
            else:
                strand = '-'
        
        if first_in_pair or not is_paired:
            if secondary or supplementary:
                if pos in generated_read.secondary:
                    generated_read.secondary[pos][2].append((pos,cigar,1,weight))
                else:
                    generated_read.secondary[pos] = [chrom,strand,[(pos,cigar,1,weight)]]
            else:
                if pos in generated_read.primary:
                    generated_read.primary[pos][2].append((pos,cigar,1,weight))
                else:
                    generated_read.primary[pos] = [chrom,strand,[(pos,cigar,1,weight)]]
        else:
            if not pair_is_mapped:
                continue
            if secondary or supplementary:
                if pairpos in generated_read.secondary:
                    generated_read.secondary[pairpos][2].append((pos,cigar,2,weight))
                else:
                    generated_read.secondary[pairpos] = [chrom,strand,[(pos,cigar,2,weight)]]
            else:
                if pairpos in generated_read.primary:
                    generated_read.primary[pairpos][2].append((pos,cigar,2,weight))
                else:
                    generated_read.primary[pairpos] = [chrom,strand,[(pos,cigar,2,weight)]]
    
    return generated_read

def read_samfile(filename, filetype=args.READTYPE, mmIO=args.mmIO):
    if mmIO:
        mmOut = {}
    infile = open(filename)
    current_ID    = None
    current_read  = None
    current_lines = []
    for line in infile:
        if line[0]=='@':
            continue
        l = line.rstrip().split('\t')
        ID = l[0]
        if ID == current_ID or current_ID is None:
            current_lines.append(line)
            current_ID = ID
        else:
            # A new read ID was encountered. Perform actions on current_read
            # and begin a new current_read with the new ID.
            if args.UNTEMP_OUT:
                current_read = generate_read_from_sam(current_lines,filetype,keep_seq=True)
            else:
                current_read = generate_read_from_sam(current_lines,filetype)
            residual = populate_bedgraphs(current_read)
            if mmIO == "out" and residual > 0:
                if residual not in mmOut:
                    mmOut[residual] = open('./mmIO.{}.{}.sam'.format(filetype,residual),'w')
                mmOut[residual].write(''.join(current_lines))
            current_lines = [line]
            current_ID = ID
    if args.UNTEMP_OUT:
        current_read = generate_read_from_sam(current_lines,filetype,keep_seq=True)
    else:
        current_read = generate_read_from_sam(current_lines,filetype)
    residual = populate_bedgraphs(current_read)
    if mmIO == "out" and residual > 0:
        if residual not in mmOut:
            mmOut[residual] = open('./mmIO.{}.{}.sam'.format(filetype,residual),'w')
        mmOut[residual].write(''.join(current_lines))
        
        for k in mmOut.keys():
            mmOut[k].close()

#############################
# CALCULATE UNIQUE COVERAGE # 
#############################

if args.mmIO == "in":
    # preload BEDGRAPH files of unique mappers
    print('Loading existing read coverage...')
    if args.SJ_OUT:
        input_name = '{}_{}_SJ.mmIO.bed'.format(args.SAMPLENAME,args.READTYPE)
        input_file = open(input_name)
        for line in input_file:
            chrom,start,end,name,count,strand = line.rstrip().split()
            count = float(count)
            i = (start,end)
            SJ[strand][chrom][i] = SJ[strand][chrom].get(i,float(0)) + count
        input_file.close()

    for STRAND in ['plus','minus']:
        if STRAND == 'plus':
            s = '+'
        elif STRAND == 'minus':
            s = '-'
        
        input_name = '{}_{}_{}.mmIO.bedgraph'.format(args.SAMPLENAME,args.READTYPE,STRAND)
        print('\t{}'.format(input_name))
        input_file = open(input_name)
        for line in input_file:
            chrom,start,end,count = line.rstrip().split()
            count = float(count)
            for i in range(int(start), int(end)):
                ENDPOINT[s][chrom][i] = ENDPOINT[s][chrom].get(i,float(0)) + count
        input_file.close()
        
        input_name = '{}_{}_{}_coverage.mmIO.bedgraph'.format(args.SAMPLENAME,args.READTYPE,STRAND)
        print('\t{}'.format(input_name))
        input_file = open(input_name)
        for line in input_file:
            chrom,start,end,count = line.rstrip().split()
            count = float(count)
            for i in range(int(start), int(end)):
                COVERAGE[s][chrom][i] = COVERAGE[s][chrom].get(i,float(0)) + count
        input_file.close()
        
        if args.UNTEMP_OUT:
            input_name = '{}_{}_{}_untemp.mmIO.bedgraph'.format(args.SAMPLENAME,args.READTYPE,STRAND)
            print('\t{}'.format(input_name))
            input_file = open(input_name)
            multikey = input_file.readline().rstrip().split('\t')[3:]
            for line in input_file:
                if line[0] == '#':
                    continue
                l = line.rstrip().split()
                chrom = l[0]
                start = int(l[1])
                end = int(l[2])
                countvals = [float(i) for i in l[3:]]
                for pos in range(start, end):
                    for i in range(len(multikey)):
                        if countvals[i] > 0:
                            if pos not in UNTEMP[s][chrom]:
                                UNTEMP[s][chrom][pos] = {}
                            UNTEMP[s][chrom][pos][multikey[i]] = UNTEMP[s][chrom][pos].get(multikey[i],0) + countvals[i]
            
            input_file.close()
    
    mmIO_files = [i for i in os.listdir('.') if i.startswith('mmIO.{}'.format(args.READTYPE))]
    args.SAMFILES = [v for k,v in sorted([(int(k.split('.')[-2]),k) for k in mmIO_files])]
    print('{} mmIO SAM files found.'.format(len(args.SAMFILES)))


if args.WRITE_BED:
    temp_bed = open('{}.{}.temp.bed'.format(args.SAMPLENAME, args.READTYPE,), 'w')

for samfile in args.SAMFILES:
    read_samfile(samfile,args.READTYPE)
    if args.mmIO == "in":
        mmnums = dict(
            [
                (
                    k,
                    sum(
                        [
                            b[1] for b in flatten(
                                [a.values() for a in v.values()]
                            )
                        ]
                    )
                ) for k,v in MM.items()
            ]
        )
        mmreadnum = sum(mmnums.values())
        multiplicity = int(samfile.split('.')[-2])
        print('Assigning {} ({} reads)'.format(samfile,mmreadnum))
        residuals = run_mm_assignment(multiplicity)
        print('    {} residuals'.format(residuals))


# Begin to populate BEDGRAPH dicts with non-unique reads using the algorithm
# defined in the function assign_multimapper().
# Within a multiplicity group: 
# 1) By start position (oscillate left->right / right->left),
#    calculate existing coverage within the mapped region of
#    reads of that multiplicity. If there are existing stranded
#    reads in at least one location, assign fractional reads
#    proportional to existing coverage at each position.
#    Otherwise, if 'allow_nonstranded', then the opposite strand
#    is also checked for coverage at this step.
# 2) Temporarily avoid naive assignment if no coverage exists.
#    Repeat the oscillation process until no more weighted
#    assignments can be made in that multiplicity group.

mmnums = dict(
    [
        (
            k,
            sum(
                [
                    b[1] for b in flatten(
                        [a.values() for a in v.values()]
                    )
                ]
            )
        ) for k,v in MM.items()
    ]
)

mmreadnum = sum(mmnums.values())
if mmreadnum > 0:
    print('\nAssigning {} multimappers: {}'.format(mmreadnum,args.READTYPE))
    residuals = 0
    for multiplicity in sorted(list(MM.keys())):
        residuals += run_mm_assignment(multiplicity)
    
    print(
        '\tRescued {}% of multimappers.'.format(
            round(
                float(sum(mmnums.values()) - residuals)/sum(mmnums.values()),
                3
            )*100
        )
    )
else:
    if args.mmIO == "out":
        print('Multimappers written to SAM files.')
    else:
        print('No residual multimappers.')

print("All reads assigned. Writing files...")

if args.SJ_OUT:
    if args.mmIO == "out":
        output_filename='{}_{}_SJ.mmIO.bed'.format(args.SAMPLENAME,args.READTYPE)
    else:
        output_filename='{}_{}_SJ.bed'.format(args.SAMPLENAME,args.READTYPE)
    
    outfile = open(output_filename, 'w')
    name = '.'
    for strand in ['+','-','.']:
        for chrom in SJ[strand].keys():
            for key in SJ[strand][chrom].keys():
                left,right = key
                count = round(SJ[strand][chrom][key],2)
                outfile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    chrom,
                    left,
                    right,
                    name,
                    count,
                    strand
                ))
    outfile.close()

for strand in ['+','-']:
    print('\t{} {}'.format(args.READTYPE,strand))
    if strand == '+':
        STRAND = 'plus'
    if strand == '-':
        STRAND = 'minus'
    #input, output_filename, digits=2, parallel=False, multi_key=False)
    if args.mmIO == "out":
        write_bedgraph_from_dict(
            ENDPOINT[strand],
            output_filename='{}_{}_{}.mmIO.bedgraph'.format(args.SAMPLENAME,args.READTYPE,STRAND),
            parallel=args.PARALLEL,
            multi_key=False
        )
        write_bedgraph_from_dict(
            COVERAGE[strand],
            output_filename='{}_{}_{}_coverage.mmIO.bedgraph'.format(args.SAMPLENAME,args.READTYPE,STRAND),
            parallel=args.PARALLEL,
            multi_key=False
        )
        if args.SIZE_CLASSES:
            write_bedgraph_from_dict(
                READLENGTHS[strand],
                output_filename='{}_{}_{}_sizeclass.mmIO.bedgraph'.format(args.SAMPLENAME,args.READTYPE,STRAND),
                parallel=args.PARALLEL,
                multi_key=args.SIZE_CLASSES
            )
        if args.UNTEMP_OUT:
            write_bedgraph_from_dict(
                UNTEMP[strand],
                output_filename='{}_{}_{}_untemp.mmIO.bedgraph'.format(args.SAMPLENAME,args.READTYPE,STRAND),
                parallel=args.PARALLEL,
                multi_key=args.UNTEMP_OUT
            )
    else:
        write_bedgraph_from_dict(
            ENDPOINT[strand],
            output_filename='{}_{}_{}.bedgraph'.format(args.SAMPLENAME,args.READTYPE,STRAND),
            parallel=args.PARALLEL,
            multi_key=False
        )
        write_bedgraph_from_dict(
            COVERAGE[strand],
            output_filename='{}_{}_{}_coverage.bedgraph'.format(args.SAMPLENAME,args.READTYPE,STRAND),
            parallel=args.PARALLEL,
            multi_key=False
        )
        if args.SIZE_CLASSES:
            write_bedgraph_from_dict(
                READLENGTHS[strand],
                output_filename='{}_{}_{}_sizeclass.bedgraph'.format(args.SAMPLENAME,args.READTYPE,STRAND),
                parallel=args.PARALLEL,
                multi_key=args.SIZE_CLASSES
            )
        if args.UNTEMP_OUT:
            write_bedgraph_from_dict(
                UNTEMP[strand],
                output_filename='{}_{}_{}_untemp.bedgraph'.format(args.SAMPLENAME,args.READTYPE,STRAND),
                parallel=args.PARALLEL,
                multi_key=args.UNTEMP_OUT
            )


print("Coverage calculations complete!")



