import sys
import os
import re
import copy
import argparse
import math
import gff_utils as gu
import fasta_utils as fu
import bedgraph_utils as bu
from collections import Counter

########################
### ARGUMENT PARSING ###
########################

parser = argparse.ArgumentParser()
parser.add_argument('-A','--annotation',dest='ANNOTATION',
                    help="Path to a reference annotation GTF/GFF.")
parser.add_argument('-N',dest='NONSTRANDED_BEDGRAPH',
                    help="Path to nonstranded bedgraph file(s).",
                    default=[], nargs = '+')
parser.add_argument('-P',dest='PLUS_BEDGRAPH',
                    help="Path to nonstranded bedgraph file(s).",
                    default=[], nargs = '+')
parser.add_argument('-M',dest='MINUS_BEDGRAPH',
                    help="Path to nonstranded bedgraph file(s).",
                    default=[], nargs = '+')
parser.add_argument("--gene", dest='GENE',
                    help="Output quantification on the gene level.",
                    default=False, action='store_true')
parser.add_argument("--buffer", dest='BUFFER',
                    help="Number of allowed nonoverlapping terminal \
                    nucleotides.", default=0, type=int)
parser.add_argument("--min_prior", dest='MIN_PRIOR',
                    help="Minimum nucleotide length to count \
                    as a prior.", default=20, type=int)
parser.add_argument("--norm", dest='NORM',
                    help="Number of allowed nonoverlapping terminal \
                    nucleotides.", default='TPM', type=str,
                    choices=['TPM','CPM','counts'])
args = parser.parse_args()

#################
### FUNCTIONS ###
#################

def which(x,value=True):
    """Returns a list of locations in x that satisfy value"""
    return [a for a,b in enumerate(x) if b==value]


def notwhich(x,value=0):
    """Returns a list of locations in x that do not satisty value"""
    return [a for a,b in enumerate(x) if b!=value]


def flatten(list_of_lists):
    """Collapses a list/tuple of lists into a single list"""
    return [item for sublist in list_of_lists for item in sublist]

def std(list):
    """Calculates the standard deviation of a list of numbers"""
    if len(list) <= 1:
        return 0
    
    list = [float(i) for i in list]
    m = sum(list)/len(list)
    diffsq = [(i - m)**2 for i in list]
    return math.sqrt(sum(diffsq)/(len(list) - 1))

def is_inside(range_a,range_b):
    """Determines if one start/end double is contained
    within the range of the second"""
    return range_a[0] >= range_b[0] and range_a[1] <= range_b[1]

def convert_to_gene(dictionary, sep='.'):
    """ Collapse the values in a dictionary
    based on the key's prefix """
    keys = dictionary.keys()
    outdict = {}
    prefixes = list(set([k.split(sep)[0] for k in keys]))
    for i in range(len(prefixes)):
        combine = [k for k in keys if prefixes[i] in k]
        cval = 0
        for k in combine:
            cval += dictionary[k]
        outdict[prefixes[i]] = cval
    return outdict

def update_priors(edge_value,edge_length,assign_to):
    """ Updates the dicts 'priors' and 'p_lens'
    for a transcript during quantification 
    with the values of a new edge """
    if len(assign_to) > 1:
        # Assign partial values proportionally
        # to each plus-stranded transcript
        expected_reads = []
        for t in assign_to:
            if t in priors:
                expected_reads += [priors[t]*edge_length]
            else:
                expected_reads += [None]
        
        sum_p = sum([i for i in expected_reads if i])
        if sum_p >= edge_value:
            # Assign the whole edge value to existing priors
            scale = float(edge_value)/sum_p
            for i in range(len(assign_to)):
                if expected_reads[i]:
                    t = assign_to[i]
                    # Update priors with a length-weighted average of the new edge
                    priors[t] = ((priors[t]*p_lens[t]) + (expected_reads[i]*scale)) / (p_lens[t] + edge_length)
                    # Add the edge length to p_lens
                    p_lens[t] += edge_length
        else:
            number_naive = len([i for i in expected_reads if i is None])
            if number_naive > 0:
                # Assign the expected weights to existing priors
                # and the residual evenly to the edges without priors
                for i in range(len(assign_to)):
                    residuals = edge_value
                    if expected_reads[i]:
                        t = assign_to[i]
                        # Update priors with a length-weighted average of the new edge
                        priors[t] = ((priors[t]*p_lens[t]) + expected_reads[i]) / (p_lens[t] + edge_length)
                        # Add the edge length to p_lens
                        p_lens[t] += edge_length
                        residuals = residuals - expected_reads[i]
                
                residual_to_assign = residuals / number_naive
                # Loop through a final time and assign residuals
                for i in range(len(assign_to)):
                    if expected_reads[i] is None:
                        t = assign_to[i]
                        priors[t] = residual_to_assign / edge_length
                        # Add the edge length to p_lens
                        p_lens[t] = edge_length
            else:
                # No naive transcripts exist, so assign the entire edge weight
                # proportionally to the existing edges
                if sum_p > 0:
                    proportions = [float(i)/sum_p for i in expected_reads]
                else:
                    proportions = [float(1)/len(assign_to) for i in expected_reads]
                
                for i in range(len(assign_to)):
                    t = assign_to[i]
                    # Update priors with a length-weighted average of the new edge
                    priors[t] = ((priors[t]*p_lens[t]) + (proportions[i]*edge_value)) / (p_lens[t] + edge_length)
                    # Add the edge length to p_lens
                    p_lens[t] += edge_length
    
    else:
        # Assign the full edge_value to the transcript
        t = assign_to[0]
        priors[t] = ((priors.get(t,0)*p_lens.get(t,0)) + edge_value) / (p_lens.get(t,0) + edge_length)
        p_lens[t] = p_lens.get(t,0) + edge_length


#############################
### IMPORT REFERENCE DATA ###
#############################

# Initialize a dictionary to store coverage data
# from all input bedgraph files
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

total_readcount = round(
    sum(
        [sum(
            [sum(i.values()) for i in coverage[strand].values()]
        ) for strand in ['+','-','.']
        ]
    )
)
    
transcripts = gu.parse_annotation(args.ANNOTATION)
IDs = sorted(list(transcripts.keys()))
print('# {} reference transcripts: {}'.format(
    len(IDs),
    args.ANNOTATION
))

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

print('# {} nonoverlapping loci'.format(
    sum([len(overlap_groups[i]) for i in overlap_groups.keys()])
))

#################################
### EVALUATE THE GIVEN SAMPLE ###
#################################

# Populate a dictionary with quantification
# values for each ID
values = {}
for chrom in overlap_groups.keys():
    for group in overlap_groups[chrom]:
        # Get the start/stop position of all
        # exons in the overlapping group
        node_dict = {}
        for ID in group:
            # Make node pairs that represent all edges
            # of a transcript model
            starts = [i - 1 for i in transcripts[ID].get_exon_start()]
            ends = transcripts[ID].get_exon_end()
            if args.BUFFER:
                # Add buffer zones as edges connected to
                # the terminal exons of the model
                buffer_left = (starts[0] - args.BUFFER, starts[0])
                buffer_right = (ends[-1], ends[-1] + args.BUFFER)
                node_dict[ID] = [buffer_left] + list(zip(starts,ends)) + [buffer_right]
            else:
                node_dict[ID] = list(zip(starts,ends))

        # An array of all unique nodes in the group
        all_nodes = sorted(list(set(flatten(flatten(node_dict.values())))))
        # Make an array of edges, represented by a double of (start,end)
        all_edges = list(zip(all_nodes[:-1],all_nodes[1:]))
        t_lens = dict([(k,sum([b-a for a,b in v])) for k,v in node_dict.items()])

        # Determine which IDs occupy each edge
        edge_assignment = {}
        for e in range(len(all_edges)):
            edge_assignment[e] = []
            for ID in node_dict.keys():
                # Determine which transcript models each edge is compatible with
                if any([is_inside(all_edges[e],n_e) for n_e in node_dict[ID]]):
                    edge_assignment[e].append(ID)
        
        # Calculate the number of reads mapping to each edge
        # for the +, -, and . data
        edge_weight_p = {}
        edge_weight_m = {}
        edge_weight_n = {}
        for e in range(len(all_edges)):
            # If the edge belongs to at least one ID,
            # get the summed coverage values for that edge
            if edge_assignment[e]:
                if chrom in coverage['+']:
                    pcov = sum([
                        coverage['+'][chrom].get(i,0)
                        for i in range(all_edges[e][0],all_edges[e][1])
                    ])
                    edge_weight_p[e] = pcov
                
                if chrom in coverage['-']:
                    mcov = sum([
                        coverage['-'][chrom].get(i,0)
                        for i in range(all_edges[e][0],all_edges[e][1])
                    ])
                    edge_weight_m[e] = mcov
                
                if chrom in coverage['.']:
                    ncov = sum([
                        coverage['.'][chrom].get(i,0)
                        for i in range(all_edges[e][0],all_edges[e][1])
                    ])
                    edge_weight_n[e] = ncov

        # Calculate priors based on unique edge weights
        unique_edges = [k for k,v in edge_assignment.items() if len(v) == 1]
        priors = {}
        p_lens = {}

        for e in unique_edges:
            # Get the unique transcript
            t = edge_assignment[e][0]
            p_lens[t] = p_lens.get(t,0) + (all_edges[e][1] - all_edges[e][0])
            # Update prior for the transcript with any nonstranded data
            priors[t] = priors.get(t,0) + edge_weight_n.get(e,0)
            # Update prior for the transcript with stranded data
            if transcripts[t].strand == '+':
                priors[t] = priors.get(t,0) + edge_weight_p.get(e,0)
            elif transcripts[t].strand == '-':
                priors[t] = priors.get(t,0) + edge_weight_m.get(e,0)
        
        # Length-normalize priors
        tdel = []
        for t in priors.keys():
            if p_lens[t] < args.MIN_PRIOR:
                tdel += [t]
            else:
                priors[t] = float(priors[t])/p_lens[t]
        for t in tdel:
            del priors[t]
            del p_lens[t]
        
        # Update priors proportionally to the unique coverage
        unassigned_edges = [k for k,v in edge_assignment.items() if len(v) > 1]
        if not priors:
            allow_naive = True
        else:
            allow_naive = False
        
        while unassigned_edges:
            trash = []
            for i in range(len(unassigned_edges)):
                e = unassigned_edges[i]
                # If there is <= 1 unknown on the edge, resolve
                # the edge using existing priors. Otherwise
                if allow_naive or sum([t not in priors for t in edge_assignment[e]]) <= 1:
                    assign_edge = True
                else:
                    assign_edge = False
                
                if assign_edge:
                    e_trans = edge_assignment[e]
                    e_len = all_edges[e][1] - all_edges[e][0]
                    e_strand = [transcripts[t].strand for t in e_trans]
                    plus_transcripts = [t for t,s in zip(e_trans,e_strand) if s == '+']
                    minus_transcripts = [t for t,s in zip(e_trans,e_strand) if s == '-']
                    # Assign plus-stranded coverage
                    if plus_transcripts:
                        e_val = edge_weight_p.get(e,0)
                        if e_val:
                            update_priors(e_val,e_len,plus_transcripts)
                    
                    # Assign minus-stranded coverage
                    if minus_transcripts:
                        e_val = edge_weight_m.get(e,0)
                        if e_val:
                            update_priors(e_val,e_len,minus_transcripts)
                    
                    # Assign nonstranded coverage proportionally
                    # based on updated priors
                    e_val = edge_weight_n.get(e,0)
                    if e_val:
                        update_priors(e_val,e_len,e_trans)
                    
                    trash += [unassigned_edges[i]]
            # Remove all the edges that were assigned
            unassigned_edges = [i for i in unassigned_edges if i not in trash]
            # If no edges were assigned in the first cycle,
            # allow naively-weighted assignment
            if len(trash) == 0:
                allow_naive = True
        
        # Length-normalize priors if the output will be in TPM
        if args.NORM in ['CPM','counts']:
            for t in priors.keys():
                priors[t] = priors[t] * t_lens[t]
        
        if args.GENE:
            priors = convert_to_gene(priors)
        
        values.update(priors)
        # for k,v in priors.items():
            # print('{}\t{}'.format(k,v))


if args.GENE:
    # Group all transcripts with the same prefix
    output_IDs = sorted(list(set([i.split('.')[0] for i in IDs])))
else:
    output_IDs = IDs
            
if args.NORM == 'TPM':
    vsum = sum(values.values())
    
for ID in output_IDs:
    v = values.get(ID, 0)
    if args.NORM == 'CPM':
        v = round(v * 10**6 / total_readcount,2)
    elif args.NORM == 'TPM':
        v = round(v * 10**6 / vsum,2)
    else:
        v = round(v)
    print('{}\t{}'.format(ID,v))
