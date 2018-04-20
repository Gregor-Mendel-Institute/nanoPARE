import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('INFILE',
                    help="Path to a reference annotation GTF/GFF.")
parser.add_argument('-F','--fields',dest='FIELDS',
                    help="List of fields to check for uniqueness",
                    default=[0,1,2,5], nargs = '+')
parser.add_argument('-S','--select',dest='SELECT',
                    help="Which feature to take as representative",
                    type=str, default='upstream',
                    choices=['upstream','downstream','left','right','first','last','highscore','lowscore'])
parser.add_argument('--startline',dest='STARTLINE',
                    help="Which field (0-based) contains start position",
                    type=int, default=1)
parser.add_argument('--endline',dest='ENDLINE',
                    help="Which field (0-based) contains end position",
                    type=int, default=2)
parser.add_argument('--strandline',dest='STRANDLINE',
                    help="Which field (0-based) contains strand",
                    type=int, default=5)
parser.add_argument('--scoreline',dest='SCORELINE',
                    help="Which field (0-based) contains score",
                    type=int, default=4)
parser.add_argument('--only_unique',dest='ONLY_UNIQUE',
                    help="Suppress writing of all features considered duplicates",
                    action='store_true', default=False)

args = parser.parse_args()
args.FIELDS = [int(i) for i in args.FIELDS]

def get_unique(array):
    ''' Takes an array of values
    and returns a unique array in
    the original order '''
    prev = {}
    unique = []
    for item in array:
        if item in prev:
            continue
        prev[item] = True
        unique.append(item)
    return unique


def output_best(bed_lines,select,startline,endline,strandline,scoreline):
    ''' Taks an array of BED file lines
    and outputs either the first or the last
    line depending on the strand '''
    if len(bed_lines) == 1:
        return bed_lines[0]
    else:
        if select == 'first':
            return bed_lines[0]
        elif select == 'last':
            return bed_lines[-1]
        
        if 'score' in select:
            # Sort based on the score value
            score = [float(i.split('\t')[scoreline]) for i in bed_lines]
            sorted_order = [b for a,b in sorted([(b,a) for a,b in enumerate(score)])]
            if select == 'highscore':
                return bed_lines[sorted_order[-1]]
            elif select == 'lowscore':
                return bed_lines[sorted_order[0]]
        
        # Sort lines by leftmost (sub-sort rightmost) position
        starts = [int(i.split('\t')[startline]) for i in bed_lines]
        ends = [int(i.split('\t')[endline]) for i in bed_lines]
        startend = list(zip(starts,ends))
        sorted_order = [b for a,b in sorted([(b,a) for a,b in enumerate(startend)])]
        
        if select == 'left':
            return bed_lines[sorted_order[0]]
        elif select == 'right':
            return bed_lines[sorted_order[-1]]
        
        strand = bed_lines[0].split('\t')[strandline]
        if select == 'upstream':
            if strand == '+':
                return bed_lines[sorted_order[0]]
            elif strand == '-':
                return bed_lines[sorted_order[-1]]
        elif select == 'downstream':
            if strand == '+':
                return bed_lines[sorted_order[-1]]
            elif strand == '-':
                return bed_lines[sorted_order[0]]
            


last_line = None
duplicate_lines = []
for line in open(args.INFILE):
    line = line.rstrip()
    # Make a hash of all selected fields for the current line
    new_feat = '_'.join([line.split('\t')[i] for i in args.FIELDS])
    if last_line:
        # Check if the new hash is the same as the last hash
        if new_feat == last_feat:
            # Store both the old and new lines in duplicate_lines
            duplicate_lines += [last_line]
            duplicate_lines += [line]
            continue
        else:
            # Output a single line based on the strand
            if not args.ONLY_UNIQUE or (args.ONLY_UNIQUE and len(duplicate_lines) == 1):
                print(output_best(get_unique(duplicate_lines),args.SELECT,args.STARTLINE,args.ENDLINE,args.STRANDLINE,args.SCORELINE))
    
    duplicate_lines = [line]
    last_feat = new_feat
    last_line = line

if not args.ONLY_UNIQUE or (args.ONLY_UNIQUE and len(duplicate_lines) == 1):
    print(output_best(get_unique(duplicate_lines),args.SELECT,args.STARTLINE,args.ENDLINE,args.STRANDLINE,args.SCORELINE))
