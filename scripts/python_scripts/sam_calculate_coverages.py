import os
import re
import sys
import argparse
import time
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
    choices=['BODY','TSS','TES'], default="BODY", type=str
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
    "--map_softclip", dest='MAP_SOFTCLIP',
    help="If true, maps 5' end regardless of softclipping",
    default=False, action="store_true"
)
parser.add_argument(
    "--minmatch", dest='MINMATCH',
    help="Minimum number of matching nucleotides for a read.",
    default=15, type=int
)

args = parser.parse_args()


######################################
# ENVIRONMENT SETUP: DATA STRUCTURES #
######################################

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

file.close()
genome_size = sum(chromosomes.values())

# The dicts COVERAGE and ENDPOINT will store read count information
# to be written as BEDGRAPH files.
print('Constructing coverage dictionaries...')
COVERAGE = {}
ENDPOINT = {}

COVERAGE['+'] = {}
COVERAGE['-'] = {}

ENDPOINT['+'] = {}
ENDPOINT['-'] = {}

for chrom,length in chromosomes.items():
    COVERAGE['+'][chrom] = {}
    COVERAGE['-'][chrom] = {}
    ENDPOINT['+'][chrom] = {}
    ENDPOINT['-'][chrom] = {}

# One class is defined for this script:
# RNASeqRead objects contain all information for a single read / read pair
# Mappings are stored as a quad in one of the lists 'primary' or 'secondary':
#   (pos, chrom, strand, cigar)

class RNASeqRead():
    def __init__(self,ID=None,paired=False,Nmap=0,primary={},secondary={}):
        self.paired    = paired
        self.Nmap      = Nmap
        self.primary   = primary
        self.secondary = secondary
        self.ID        = ID
    
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
                positions.append(list(self.secondary.keys()))
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


def populate_bedgraphs(read_object):
    """Adds values from read_object to the appropriate bedgraph dicts."""
    if len(read_object.primary)==0:
        return
    
    # If multimapper, store in MM and do not change bedgraph dicts
    nmap = read_object.Nmap
    if nmap > 1:
        minpos = read_object.minpos()
        cigar  = read_object.mincigar()
        read_object.ID = None
        if nmap in MM:
            if minpos in MM[nmap]:
                if cigar in MM[nmap][minpos]:
                    MM[nmap][minpos][cigar] = (read_object,MM[nmap][minpos][cigar][1]+1)
                else:
                    MM[nmap][minpos][cigar] = (read_object,1)
            else:
                MM[nmap][minpos] = {cigar:(read_object,1)}
        else:
            MM[nmap] = {minpos:{cigar:(read_object,1)}}
    else:
        # Parse the mapping information of the uniquely mapping read
        chrom,strand,poslist = read_object.map_locations()[0]
        map_positions = set()
        for pos,cigar in poslist:
            mate_positions = get_mapping_positions(pos,cigar)
            mate_positions = set([i for i in mate_positions if i > 0 and i <= chromosomes[chrom]])
            map_positions.update(mate_positions)
        
        if len(map_positions) > args.MINMATCH:
            if strand == '+':
                ENDPOINT_end = min(map_positions) - 1
            elif strand == '-':
                ENDPOINT_end = max(map_positions) - 1
            
            if args.READTYPE == 'TES':
                if strand == '+':
                    strand = '-'
                else:
                    strand = '+'
            
            ENDPOINT[strand][chrom][ENDPOINT_end] = ENDPOINT[strand][chrom].get(ENDPOINT_end,0) + 1
            for pos in map_positions:
                COVERAGE[strand][chrom][pos-1] = COVERAGE[strand][chrom].get(pos-1,0) + 1

def assign_multimapper(
    read_object,
    value=1,
    allow_naive=False,
    include_secondary=args.SECONDARY,
    allow_nonstranded=args.ALLOW_NS
):
    # Assigns values to the respective bedgraphs for reads mapping to multiple locations
    # using a 'rich-get-richer' algorithm:
    # 1) If stranded data of the same read type exists within region of coverage, assign proportionally
    # 2) Else if coverage data of any type exists, assign proportionally
    # 3) Else if allow_naive, then distribute evenly across all mapping positions
    # 4) Else do nothing, returning the original read_object
    
    # Get dictionary of all mappings from the read_object
    all_mappings = read_object.map_locations()
    
    mappos_list = [
        set(
            flatten(
                [
                    [
                        i for i in list(get_mapping_positions(p,c))
                        if i > 0 and i <= chromosomes[chrom]
                    ] for p,c in poslist
                ]
            )
        )
        for chrom,strand,poslist in all_mappings
    ]
    chroms_list = [chrom for chrom,strand,poslist in all_mappings]
    long_enough = [len(i) >= args.MINMATCH for i in mappos_list]
    all_mappings = [all_mappings[i] for i in which(long_enough)]
    mappos_list = [mappos_list[i] for i in which(long_enough)]
    chroms_list = [chrom for chrom,strand,poslist in all_mappings]
    strand_list = [strand for chrom,strand,poslist in all_mappings]
    
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
    total_coverage = sum(existing_reads)
    proportions = [float(i)/total_coverage for i in existing_reads]
    
    for i in range(len(all_mappings)):
        if strand_list[i] == '+':
            ENDPOINT_end = min(mappos_list[i]) - 1
        elif strand_list[i] == '-':
            ENDPOINT_end = max(mappos_list[i]) - 1
        
        if args.READTYPE == 'TES':
            if strand_list[i] == '+':
                strand_list[i] = '-'
            else:
                strand_list[i] = '+'
        
        ENDPOINT[strand_list[i]][chroms_list[i]][ENDPOINT_end] = \
            ENDPOINT[strand_list[i]][chroms_list[i]].get(ENDPOINT_end,0) + \
            float(value)*proportions[i]
        for pos in mappos_list[i]:
            COVERAGE[strand_list[i]][chroms_list[i]][pos-1] = \
                COVERAGE[strand_list[i]][chroms_list[i]].get(pos-1,0) + \
                float(value)*proportions[i]        

def write_bedgraph_from_dict(input,output_filename,digits=2):
    """Writes unsorted BEDGRAPH to output_filename from input dict"""
    
    def generate_bedgraph_lines(values_dict,chrom,queue):
        """Converts dict of values to bedgraph format"""
        start = 0
        prevpos = 0
        prevcount = None
        chromlen = len(input[chrom])
        position = None
        # Iterate through every position with values in the dict
        for position in sorted(list(values_dict.keys())):
            count = round(values_dict[position],digits)
            if count != prevcount or int(position) > 1 + prevpos:
                # The newly encountered value is not a continuation
                # of the previous value. Write the old run and start another.
                if prevcount > 0:
                    # Write the old run to outfile
                    queue.put(
                        '\t'.join(
                            [
                                str(i) for i in [
                                    chrom,
                                    start,
                                    prevpos + 1,
                                    prevcount
                                ]
                            ]
                        ) + '\n'
                    )
                    
                start = position
            prevcount = count
            prevpos = int(position)
        
        if position and prevcount > 0:
            queue.put(
                '\t'.join(
                    [
                        str(i) for i in [
                            chrom,
                            start,
                            prevpos + 1,
                            prevcount
                        ]
                    ]
                ) + '\n'
            )
        
    queue = mp.Queue()
    STOP_TOKEN = "FINISHED"
    writer_process = mp.Process(
        target=writer,
        args=(output_filename,queue,STOP_TOKEN)
    )
    writer_process.start()
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
        elif operator in ['D','N']:
            current_pos += number
        elif operator in ['I','H']:
            continue
    
    return positions
        

def generate_read_from_sam(input_lines,readtype=args.READTYPE):
    """Convert a set of lines from a SAM file into an RNASeqRead object"""
    generated_read = RNASeqRead(ID=None,Nmap=0,primary={},secondary={})
    for line in input_lines:
        l = line.rstrip().split('\t')
        ID = l[0]
        if not generated_read.ID:
            generated_read.ID = ID
        assert ID == generated_read.ID, 'ERROR: nonmatching IDs in input_lines:\n{}\t{}'.format(generated_read.ID,ID)
        attributes = dict([(':'.join(i.split(':')[0:2]),i.split(':')[-1]) for i in l[11:]])
        Nmap       = int(attributes['NH:i'])
        if generated_read.Nmap == 0:
            generated_read.Nmap = Nmap
        assert Nmap == generated_read.Nmap, 'ERROR: inconsistent Nmap score for {}'.format(ID)
        
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
        
        generated_read.paired = is_paired
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
                    generated_read.secondary[pos][2].append((pos,cigar))
                else:
                    generated_read.secondary[pos] = [chrom,strand,[(pos,cigar)]]
            else:
                if pos in generated_read.primary:
                    generated_read.primary[pos][2].append((pos,cigar))
                else:
                    generated_read.primary[pos] = [chrom,strand,[(pos,cigar)]]
        else:
            if not pair_is_mapped:
                continue
            if secondary or supplementary:
                if pairpos in generated_read.secondary:
                    generated_read.secondary[pairpos][2].append((pos,cigar))
                else:
                    generated_read.secondary[pairpos] = [chrom,strand,[(pos,cigar)]]
            else:
                if pairpos in generated_read.primary:
                    generated_read.primary[pairpos][2].append((pos,cigar))
                else:
                    generated_read.primary[pairpos] = [chrom,strand,[(pos,cigar)]]
    
    return generated_read

def read_samfile(filename,filetype=args.READTYPE):
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
            current_read = generate_read_from_sam(current_lines,filetype)
            populate_bedgraphs(current_read)
            current_lines = [line]
            current_ID = ID
    current_read = generate_read_from_sam(current_lines,filetype)
    populate_bedgraphs(current_read)

#############################
# CALCULATE UNIQUE COVERAGE # 
#############################

for samfile in args.SAMFILES:
    print('Reading {}...'.format(samfile))
    read_samfile(samfile,args.READTYPE)


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

oscillator = False
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

print('\nAssigning {} multimappers: {}'.format(sum(mmnums.values()),args.READTYPE))
residuals = 0
for multiplicity in sorted(list(MM.keys())):
    mmnum = mmnums[multiplicity]
    print('\t{}: {}'.format(multiplicity,mmnum))
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
    residuals += sum(
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

print(
    '\tRescued {}% of multimappers.'.format(
        round(
            float(sum(mmnums.values()) - residuals)/sum(mmnums.values()),
            3
        )*100
    )
)

print("All reads assigned. Writing BEDGRAPH files...")
for strand in ['+','-']:
    print('\t{} {}'.format(args.READTYPE,strand))
    if strand == '+':
        STRAND = 'plus'
    if strand == '-':
        STRAND = 'minus'
    write_bedgraph_from_dict(
        ENDPOINT[strand],
        output_filename='{}_{}_{}.bedgraph'.format(args.SAMPLENAME,args.READTYPE,STRAND)
    )
    write_bedgraph_from_dict(
        COVERAGE[strand],
        output_filename='{}_{}_{}_coverage.bedgraph'.format(args.SAMPLENAME,args.READTYPE,STRAND)
    )

print("Coverage calculations complete!")



