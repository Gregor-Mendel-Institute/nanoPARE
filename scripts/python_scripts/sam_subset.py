import os
import re
import sys
import argparse
import time
import fasta_utils as fu
import gff_utils as gu
import multiprocessing as mp

# Takes a sorted SAM file and annotation GTF/GFF file as inputs
# Outputs only reads that overlap with the provided annotations

###################
# INPUT ARGUMENTS #
###################

parser = argparse.ArgumentParser()
parser.add_argument(
    '-A','--annotation',dest='ANNOTATION',
    help="Path to a reference annotation GTF/GFF."
)
parser.add_argument(
    "-I", "--input", dest='SAMFILES', nargs='+',
    help="Path(s) to aligned SAM file(s). Require header."
)
parser.add_argument(
    "-U", "--unique", dest='UNIQUE',
    help="If true, only outputs the first read in each unique position.",
    default=False, action="store_true"
)
parser.add_argument(
    "--secondary", dest='SECONDARY',
    help="If true, includes non-primary mappings.",
    default=False, action="store_true"
)
parser.add_argument(
    "--antisense", dest='ANTISENSE',
    help="If true, outputs reads mapping antisense to the annotations.",
    default=False, action="store_true"
)
parser.add_argument(
    "-F", "--fasta", dest='FASTA',
    help="Genome FASTA file",
    default=None, type=str
)
parser.add_argument(
    "--feature", dest='FEATURE',
    help="Feature type to plot",
    default='transcript',choices=['CDS','transcript','5UTR','3UTR'],type=str
)
parser.add_argument(
    "--buffer", dest='BUFFER',
    help="Number of nucleotides to extend each annotation",
    default=0, type=int
)
parser.add_argument(
    "--minmatch", dest='MINMATCH',
    help="Minimum number of matching nucleotides for a read.",
    default=18, type=int
)
parser.add_argument(
    "--softclip_type", dest='SOFTCLIP_TYPE',
    help="Choose which softclipped reads to keep.",
    default='both', choices=['none','5p','3p','both']
)
parser.add_argument(
    "--nucfreqs", dest='NUCFREQS',
    help="Saves a table of nucleotide frequencies in selected genome region.",
    default=False, action="store_true"
)

args = parser.parse_args()

######################################
# ENVIRONMENT SETUP: DATA STRUCTURES #
######################################

if args.FASTA:
    genome = fu.import_genome(args.FASTA)

if args.SOFTCLIP_TYPE != 'none' and not args.FASTA:
    print('ERROR: Untemplated nucleotide analysis requires a reference genome. Use the --fasta argument.')
    sys.exit(1)

chromosomes = {}

file = open(args.SAMFILES[0])
line = file.readline()
while line[0] == '@':
    print(line.rstrip())
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

# The dicts COVERAGE store
ANNO = {}
ANNO['+'] = {}
ANNO['-'] = {}

COVERAGE = {}
COVERAGE['+'] = {}
COVERAGE['-'] = {}
for chrom,length in chromosomes.items():
    ANNO['+'][chrom] = set()
    ANNO['-'][chrom] = set()
    COVERAGE['+'][chrom] = set()
    COVERAGE['-'][chrom] = set()


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

################################
# ENVIRONMENT SETUP: FUNCTIONS #
################################
    
def flatten(list_of_lists):
    return [item for sublist in list_of_lists for item in sublist]


def which(x,value=True):
    return [a for a,b in enumerate(x) if b==value]


def get_mapping_positions(pos,cigar,map_softclip=False):
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
            if map_softclip:
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
            # TODO: Add donor and acceptor sites to read_object mapping
            leftside = current_pos
            current_pos += number
            rightside = current_pos - 1
            junctions += [(leftside,rightside)]
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


def generate_read_from_sam(input_lines,keep_seq=True):
    """Convert a set of lines from a SAM file into an RNASeqRead object"""
    generated_read = RNASeqRead(ID=None,Nmap=0,primary={},secondary={},seq1=None,seq2=None)
    for line in input_lines:
        l = line.rstrip().split('\t')
        ID = l[0]
        if not generated_read.ID:
            generated_read.ID = ID
        assert ID == generated_read.ID, 'ERROR: nonmatching IDs in input_lines:\n{}\t{}'.format(generated_read.ID,ID)
        
        attributes = dict([(':'.join(i.split(':')[0:2]),i.split(':')[-1]) for i in l[11:]])
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
                    generated_read.secondary[pos][2].append((pos,cigar,1))
                else:
                    generated_read.secondary[pos] = [chrom,strand,[(pos,cigar,1)]]
            else:
                if pos in generated_read.primary:
                    generated_read.primary[pos][2].append((pos,cigar,1))
                else:
                    generated_read.primary[pos] = [chrom,strand,[(pos,cigar,1)]]
        else:
            if not pair_is_mapped:
                continue
            if secondary or supplementary:
                if pairpos in generated_read.secondary:
                    generated_read.secondary[pairpos][2].append((pos,cigar,2))
                else:
                    generated_read.secondary[pairpos] = [chrom,strand,[(pos,cigar,2)]]
            else:
                if pairpos in generated_read.primary:
                    generated_read.primary[pairpos][2].append((pos,cigar,2))
                else:
                    generated_read.primary[pairpos] = [chrom,strand,[(pos,cigar,2)]]
    
    return generated_read


def check_if_in_subset(chrom,strand,endpos):
    """ Checks whether the 5P position of a read
    is contained in the provided annotations. """
    if endpos in ANNO[strand][chrom]:
        return True
    else:
        return False

def check_if_unique(chrom,strand,endpos):
    """ Checks whether the 5P position of a read
    is already represented in the COVERAGE dict.
    Returns boolean. """
    if endpos not in COVERAGE[strand][chrom]:
        COVERAGE[strand][chrom].add(endpos)
        return True
    else:
        return False


def read_samfile(filename):
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
            if args.SOFTCLIP_TYPE == 'none':
                current_read = generate_read_from_sam(current_lines)
            else:
                current_read = generate_read_from_sam(current_lines,keep_seq=True)
            
            current_lines = [line]
            current_ID = ID
            if current_read.primary:
                mapping = list(current_read.primary.values())[0]
            elif args.SECONDARY and current_read.secondary:
                mapping = list(current_read.secondary.values())[0]
            else:
                continue
            chrom = mapping[0]
            strand = mapping[1]
            positions = get_mapping_positions(mapping[2][0][0],mapping[2][0][1])[0]
            if len(positions) < args.MINMATCH:
                continue
            
            if strand == '+':
                endpos = min(positions)
            else:
                endpos = max(positions)
            
            is_contained = check_if_in_subset(chrom,strand,endpos)
            if is_contained:
                if args.UNIQUE:
                    is_unique = check_if_unique(chrom,strand,endpos)
                    if is_unique:
                        for samline in current_lines:
                            print(samline.rstrip())
                else:
                    for samline in current_lines:
                        print(samline.rstrip())


############
# EVALUATE # 
############

transcripts = gu.parse_annotation(args.ANNOTATION)
IDs = sorted(list(transcripts.keys()))
new_transcripts = {}

# Construct a dictionary of leftmost->rightmost positions
# for each transcript in the annotation file
start_end_index = {}
for ID,transcript in transcripts.items():
    chrom = transcript.chrom
    start = transcript.start - args.BUFFER
    stop = transcript.stop + args.BUFFER
    if chrom not in start_end_index:
        start_end_index[chrom] = {}
    start_end_index[chrom][start] = \
        start_end_index[chrom].get(start,[]) + [(stop,ID)]

# Generate a collection of 'overlap groups' that share
# at least 1nt +-buffer
overlap_groups = {}
def get_rpos(chrom,pos):
    return max([i[0] for i in start_end_index[chrom].get(pos,None)])

def get_IDs(chrom,pos):
    return [i[1] for i in start_end_index[chrom].get(pos,None)]

for chrom in start_end_index.keys():
    # start an empty array for the chromosome
    overlap_groups[chrom] = []
    # identify the lowest and highest start positions in the chromosome
    left_positions = sorted(start_end_index[chrom].keys())
    next_loc = 1
    lpos = left_positions[0]
    last_position = left_positions[-1]
    # if only one position exists, output IDs from this position
    if last_position == lpos:
        overlap_groups[chrom].append(tuple(get_IDs(chrom,lpos)))
    else:
        # if more than one starting position exists, then
        # iterate over the ordered start positions and aggregate
        # an array of IDs reachable from that position.
        while lpos < last_position:
            # get rightmost postion of all features starting at pos
            contained_IDs = []
            contained_IDs += get_IDs(chrom,lpos)
            rpos = get_rpos(chrom,lpos)
            lnext = left_positions[next_loc]
            # as long as overlapping transcripts can be found,
            # add their IDs to the group
            while lnext <= rpos:
                contained_IDs += get_IDs(chrom,lnext)
                n_rpos = get_rpos(chrom,lnext)
                if n_rpos > rpos:
                    rpos = n_rpos
                next_loc += 1
                
                if next_loc < len(left_positions):
                    lnext = left_positions[next_loc]
                else:
                    break
            
            # when lnext is beyond rpos, dump the 
            # current group and start a new one
            overlap_groups[chrom].append(tuple(contained_IDs))
            lpos = lnext
            next_loc += 1

for chrom in overlap_groups.keys():
    for group in overlap_groups[chrom]:
        all_strands = set([transcripts[t].strand for t in group])
        if len(all_strands) > 1:
            continue
        s = list(all_strands)[0]
        if args.ANTISENSE:
            if s == '+':
                s = '-'
            elif s == '-':
                s = '+'
        
        for t in group:
            positions = set(flatten([list(range(a-1,b)) for a,b in zip(transcripts[t].get_exon_start(),transcripts[t].get_exon_end())]))
            
            if args.FEATURE != 'transcript':
                positions = sorted(list(positions))
                strand = transcripts[t].strand
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
            
            ANNO[s][transcripts[t].chrom].update(positions)

nucleotide_frequencies = {}
if args.NUCFREQS:
    for strand in ['+','-']:
        for chrom in ANNO[strand].keys():
            for nuc in ANNO[strand][chrom]:
                n = genome[chrom][nuc]
                if strand == '-':
                    n = fu.complement(n)
                nucleotide_frequencies[n] = nucleotide_frequencies.get(n,0) + 1
    total_nucs = nucleotide_frequencies['A'] + nucleotide_frequencies['C'] + nucleotide_frequencies['G'] + nucleotide_frequencies['T']
    out_table = open('nucfreqs.txt','w')
    for n in ['A','C','G','T']:
        out_table.write('{}\t{}\n'.format(n,float(nucleotide_frequencies[n])/total_nucs))
    
    out_table.close()


for samfile in args.SAMFILES:
    read_samfile(samfile)


