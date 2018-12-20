import copy
import sys
import argparse

########################
### ARGUMENT PARSING ###
########################

parser = argparse.ArgumentParser()
parser.add_argument('-I','--input',dest='INPUT',
                    help="Path(s) to input BED file(s).",
                    nargs='+',default=[])
parser.add_argument("-L",'--labels', dest='LABELS',
                    help="One classifier for each input.",
                    default=[],nargs='+')
parser.add_argument("--buffer", dest='BUFFER',
                    help="Number of allowed nonoverlapping terminal \
                    nucleotides.", default=0, type=int)
parser.add_argument("--score_column", dest='SCORE_COLUMN',
                    help="Which column (0-indexed) to use as the score.", 
                    default=4, type=int)
parser.add_argument("-S","--stranded",dest='STRANDED',
                    default=False,action='store_true')
args = parser.parse_args()

###############
### CLASSES ###
###############

class BedFeature():
    def __init__(self,chrom,start,end,strand,ID='',score=0,other=[]):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.ID = ID
        self.score = score
        self.other = other
    
    def __repr__(self):
        return '<BED feature ({}:{}-{} {})>'.format(
            self.chrom,
            self.start,
            self.end,
            self.strand
        )
    
    def overlaps(self,feature,stranded,buffer=0):
        if stranded:
            if self.strand != feature.strand:
                return False
        
        o = sum([
            feature.start - buffer >= self.start,
            feature.start - buffer >= self.end,
            feature.end + buffer >= self.start,
            feature.end + buffer >= self.end,
        ])
        
        if o == 0 or o == 4:
            return False
        
        return True
    
    def write(self):
        ''' Writes the bed feature as a
        tab-separated line '''
        line = '\t'.join([
            str(i)
            for i in [
                self.chrom,
                self.start,
                self.end,
                self.ID,
                self.score,
                self.strand,
                '\t'.join([j for j in other])
            ]
        ])
        print(line)
        


#################
### FUNCTIONS ###
#################

def which(x,value=True):
    """Returns a list of locations in x that satisfy value"""
    return [a for a,b in enumerate(x) if b==value]

def table(x):
    """Returns a list of elements and their frequency in x,
    sorted by the number of times they occur"""
    y=sorted(list(set(x)))
    z={}
    for i in y:
        z[i]=len([1 for b in x if b==i])
    return [(k,v) for v,k in sorted([(v,k) for k,v in z.items()],reverse=True)]

def notwhich(x,value=0):
    """Returns a list of locations in x that do not satisty value"""
    return [a for a,b in enumerate(x) if b!=value]


def flatten(list_of_lists):
    """Collapses a list/tuple of lists into a single list"""
    return [item for sublist in list_of_lists for item in sublist]

def is_inside(range_a,range_b):
    """Determines if one start/end double is contained
    within the range of the second"""
    return range_a[0] >= range_b[0] and range_a[1] <= range_b[1]

def get_overlap_groups(feature_dict, buffer):
    """ Identifies all groups of overlapping features.
    Returns an array of tuples."""
    # Construct a dictionary of leftmost->rightmost positions
    # for each feature in the input files
    def get_rpos(chrom,pos):
        return max([i[0] for i in start_end_index[chrom].get(pos,None)])
    
    def get_IDs(chrom,pos):
        return [i[1] for i in start_end_index[chrom].get(pos,None)]
    
    start_end_index = {}
    for ID,feature in feature_dict.items():
        chrom = feature.chrom
        start = feature.start - buffer
        end = feature.end + buffer
        if chrom not in start_end_index:
            start_end_index[chrom] = {}
        start_end_index[chrom][start] = \
            start_end_index[chrom].get(start,[]) + [(end,ID)]
    
    # Generate a collection of 'overlap groups' that share
    # at least 1nt +-buffer
    overlap_groups = {}
    
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
    
    return overlap_groups

def get_consensus_model(IDs,stranded=False):
    ''' Executes a decision tree that returns a
    feature or features that represent the consensus
    of the overlap group '''

    output_list = []
    overlap_list = []
    all_starts = [features_dict[i].start for i in IDs]
    all_ends = [features_dict[i].end for i in IDs]
    all_strands = [features_dict[i].strand for i in IDs]
    all_chroms = [features_dict[i].chrom for i in IDs]
    all_labels = [features_dict[i].label for i in IDs]
    all_scores = [features_dict[i].score for i in IDs]
    
    assert len(set(all_chroms)) == 1
    
    for a in range(len(IDs)-1):
        for b in range(a+1,len(IDs)):
            overlaps = features_dict[IDs[a]].overlaps(features_dict[IDs[b]],stranded)
            overlap_list.append(overlaps)
    
    if sum(overlap_list) == len(overlap_list):
        # If no non-overlaps exist, merge all features together
        score = sum(all_scores)
        # Get the dominant label
        label_table = table(all_labels)
        if len(label_table) == 1:
            # Only one label type exists
            label = label_table[0][0]
        else:
            if label_table[0][1] > label_table[1][1]:
                # A more common label type exists
                label = label_table[0][0]
            else:
                # Pick the label of the feature with the highest score
                label = all_labels[sorted([(v,k) for k,v in enumerate(all_scores)],reverse=True)[0][1]]
        
        if label:
            ID = '.'.join([all_chroms[0],label])
        else:
            ID = all_chroms[0]
        
        if len(set(all_strands)) == 1:
            output_feature = BedFeature(
                all_chroms[0],
                min(all_starts),
                max(all_ends),
                all_strands[0],
                ID,score
            )
        else:
            # If merging features of multiple strands, remove strand label
            output_feature = BedFeature(
                all_chroms[0],
                min(all_starts),
                max(all_ends),
                '.',
                ID,score
                )
        
        output_list.append(output_feature)
    else:
        # At least one pair of features has no overlap
        if stranded:
            # Resolve stranded subsets of the locus if applicable
            plus_features = [k for k,v in zip(IDs,all_strands) if v == '+']
            minus_features = [k for k,v in zip(IDs,all_strands) if v == '-']
            used_features = []
            if plus_features and len(plus_features) < len(IDs):
                output_list += get_consensus_model(plus_features,stranded=True)
                used_features += plus_features
            if minus_features and len(minus_features) < len(IDs):
                output_list += get_consensus_model(minus_features,stranded=True)
                used_features += minus_features
            
            IDs = [i for i in IDs if i not in used_features]
            if not IDs:
                return output_list
        
        # Resolve incomplete overlap by building an edge network
        # print("UNRESOLVED - {} - {}:{}-{}".format(
            # IDs,all_chroms[0],min(all_starts),max(all_ends)
        # ))
        start_groups = [sum([int(i >= j) for j in all_ends]) for i in all_starts]
        end_groups = [sum([int(i <= j) for j in all_starts]) for i in all_ends]
        edge_count = {}
        edge_score = {}
        start_nodes = {}
        end_nodes = {}
        for ID_num in range(len(IDs)):
            s = start_groups[ID_num]
            e = end_groups[ID_num]
            if s not in edge_score:
                edge_score[s] = {}
            edge_score[s][e] = edge_score[s].get(e,float(0)) + all_scores[ID_num]
            
            if s not in start_nodes:
                start_nodes[s] = all_starts[ID_num]
                prev_score = edge_score[s][e]
            if e not in end_nodes:
                end_nodes[e] = all_ends[ID_num]
            
            elif edge_score[s][e] > prev_score:
                start_nodes[s] = all_starts[ID_num]
                end_nodes[e] = all_ends[ID_num]
            
            if s not in edge_count:
                edge_count[s] = {}
            edge_count[s][e] = edge_count[s].get(e,0) + 1
            
        # Resolve the network from left to right based first on edge count,
        # and secondarily on edge score in the case of a tie.
        furthest = 0
        for i in sorted(list(edge_count.keys())):
            if start_nodes[i] <= furthest:
                # Don't allow overlap with previously chosen feature
                del edge_count[i]
                del edge_score[i]
                # print('# removed internal start at {}:{}'.format(all_chroms[0],start_nodes[i]))
                continue
            
            choices = sorted([(v,k) for k,v in edge_count[i].items()],reverse=True)
            if choices:
                # Edges exist out of this node
                if len(choices) == 1:
                    # Only one edge exists, pick it
                    picked = choices[0][1]
                else:
                    # Multiple edges exist
                    if choices[0][0] > choices[1][0]:
                        # One dominant edge exists, pick it
                        picked = choices[0][1]
                    else:
                        # Edge counts are tied, check downstream
                        tied_ends = [k for v,k in choices if v == choices[0][0]]
                        downstream_edges = flatten([
                            [
                                (k,v)
                                for v in edge_count[k].keys()
                                if v <= max(tied_ends)
                                and edge_count[k][v] > choices[0][0]
                            ]
                            for k in edge_count.keys() if k > i
                        ])
                        if downstream_edges:
                            # Store the leftmost upstream end and
                            # pick the dominant downstream edge
                            picked = sorted([(end_nodes[k],k) for k in tied_ends])[0][1]
                            
                            if all_labels:
                                # Get the dominant label
                                sublabels = [
                                    all_labels[l] for l in range(len(all_labels))
                                    if start_groups[l] == i and end_groups[l] == picked
                                ]
                                subscores = [
                                    all_scores[l] for l in range(len(all_scores))
                                    if start_groups[l] == i and end_groups[l] == picked
                                ]
                                label_table = table(sublabels)
                                if len(label_table) == 1:
                                    # Only one label type exists
                                    label = label_table[0][0]
                                else:
                                    if label_table[0][1] > label_table[1][1]:
                                        # A more common label type exists
                                        label = label_table[0][0]
                                    else:
                                        # Pick the label of the feature with the highest score
                                        label = sublabels[sorted([(v,k) for k,v in enumerate(subscores)],reverse=True)[0][1]]
                                
                                ID = '.'.join([all_chroms[0],label])
                            else:
                                ID = all_chroms[0]
                            
                            output_feature = BedFeature(
                                all_chroms[0],
                                start_nodes[i],
                                end_nodes[picked],
                                all_strands[ID_num],
                                ID,edge_score[i][picked]
                            )
                            output_list.append(output_feature)
                            # Remove the starting edge from the network
                            del edge_count[i]
                            del edge_score[i]
                            # Pick the first dominant downstream edge
                            i = downstream_edges[0][0]
                            picked = downstream_edges[0][1]
                            furthest = end_nodes[picked]
                        else:
                            # Pick based on score
                            choices = sorted([(edge_score[i][e],e) for e in tied_ends],reverse=True)
                            picked = choices[0][1]
                            furthest = end_nodes[picked]
            
            if all_labels:
                # Get the dominant label for the chosen edge
                sublabels = [
                    all_labels[l] for l in range(len(all_labels))
                    if start_groups[l] == i and end_groups[l] == picked
                ]
                subscores = [
                    all_scores[l] for l in range(len(all_scores))
                    if start_groups[l] == i and end_groups[l] == picked
                ]
                label_table = table(sublabels)
                if len(label_table) == 1:
                    # Only one label type exists
                    label = label_table[0][0]
                else:
                    if label_table[0][1] > label_table[1][1]:
                        # A more common label type exists
                        label = label_table[0][0]
                    else:
                        # Pick the label of the feature with the highest score
                        label = sublabels[sorted([(v,k) for k,v in enumerate(subscores)],reverse=True)[0][1]]
                
                ID = '.'.join([all_chroms[0],label])
            else:
                ID = all_chroms[0]
            
            output_feature = BedFeature(
                all_chroms[0],
                start_nodes[i],
                end_nodes[picked],
                all_strands[ID_num],
                ID,edge_score[i][picked]
            )
            output_list.append(output_feature)
            furthest = end_nodes[picked]
        
    return output_list
    

#############################
### IMPORT REFERENCE DATA ###
#############################

# Initialize a dictionary to store feature data
# from all input bed files
features_dict = {}
running_tag = 0
label = None
for i in range(len(args.INPUT)):
    file = open(args.INPUT[i])
    if args.LABELS:
        label = args.LABELS[i]
    
    for line in file:
        if line[0] == '#':
            continue
        l = line.rstrip().split('\t')
        chrom = l[0]
        start = int(l[1])
        end = int(l[2])
        ID = l[3]
        score = float(l[4])
        strand = l[5]
        if len(l) > 6:
            other = l[6:]
        else:
            other = []
        
        feature = BedFeature(
            chrom,start,end,strand,ID,score,other
        )
        feature.label = label
        
        if running_tag in features_dict:
            print("ERROR: tag already exists")
            sys.exit(1)
        
        features_dict[running_tag] = feature
        running_tag += 1


overlap_groups = get_overlap_groups(features_dict, args.BUFFER)

print('# {} nonoverlapping loci'.format(
    sum([len(overlap_groups[i]) for i in overlap_groups.keys()])
))

#################################
### EVALUATE THE GIVEN SAMPLE ###
#################################

# Execute the consensus split/merge decisions
# on each overlap group

values = {}

for chrom in sorted(list(overlap_groups.keys())):
    running_tag = 0
    for overlap_IDs in overlap_groups[chrom]:
        consensus = get_consensus_model(overlap_IDs,stranded=args.STRANDED)
        if consensus:
            for feature in consensus:
                feature.ID = feature.ID + '.' + str(running_tag)
                feature.write()
                running_tag += 1
