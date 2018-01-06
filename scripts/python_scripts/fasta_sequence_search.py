############################
# GLOBAL ENVIRONMENT SETUP #
############################

import sys
import re
import os
import time
import argparse
import fasta_utils as fu
import multiprocessing as mp

########################
### ARGUMENT PARSING ###
########################

parser = argparse.ArgumentParser()
parser.add_argument(dest='genome_FASTA',
                    help="Path to a FASTA file.")
parser.add_argument(dest='match_file',
                    help="Path to a table of values to search for (see description).")
parser.add_argument('-O', '--out_folder', type=str, default='./',
                    help="Folder to output resulting BED files.")
parser.add_argument('--parallel', action='store_true',
                    help="Multiprocessing of each element in genome_FASTA")
parser.add_argument("--sort_results", action='store_true',
                    help="Whether to sort results prior to writing.")

args = parser.parse_args()

out_folder = args.out_folder
if out_folder[-1] != '/':
    out_folder+='/'

def flatten(list_of_lists):
    return [item for sublist in list_of_lists for item in sublist]

def distance(sequencelist,dist,wildcard=False):
    if type(sequencelist) is str:
        sequencelist = [sequencelist]
    for i in range(dist):
        if wildcard:
            sequencelist = list(set(flatten(
                [
                    [
                    '.'.join([string[0:j],string[j+1:]])
                    for j in range(len(string))
                    ]
                for string in sequencelist
                ]
            )))
        else:
            sequencelist = list(set(flatten(
                [
                    [
                        (IUPACnot[string[j]]).join([string[0:j],string[j+1:]])
                        for j in range(len(string))
                    ]
                for string in sequencelist
                ]
            )))
    return sequencelist


def lookahead(string):
    """Adds a lookahead assertion to the first character of a regex."""
    if string[0] != '[':
        cutposition = 0
    else:
        cutposition = min([i for i in range(len(string)) if string[i]==']'])
    return string[0:cutposition+1] + '(?=' + string[cutposition+1:] + ')'

def rc(sequence):
    return ''.join(reversed([IUPACcomp[i] for i in sequence.upper()]))

def to_regex(sequence):
    return ''.join([IUPAChash[i] for i in sequence.upper()])

def writer(merge_filename,shared_queue,stop_token):
    """Initializes a writer to output the mp queue."""
    dest_file = open(merge_filename,'w')
    while True:
        line = shared_queue.get()
        if line == stop_token:
            dest_file.close()
            return
        dest_file.write(line)

def IUPAC_matches(sequence,dist,chrom,queue,parallel):
    """Passes BED-formatted lines of all matches to an IUPAC sequence to queue.
    
    Usage:
        sequence = String following IUPAC nucelotide nomenclature
        dist = Number of allowed mismatches
        chrom = Key of chromosome in genome dict
        queue = Which queue to add BED lines
    """
    print('{} ({}) started.'.format(chrom,sequence))
    sequence = sequence.split(',')
    seqlen   = [len(i) for i in sequence]
    if type(dist) is str:
        dist = [int(i) for i in dist.split(',')]
    elif type(dist) is int:
        dist = [dist]*len(sequence)
    search_sequences = []
    assert len(sequence) == len(dist)
    for s,d in zip(sequence,dist):
        for i in range(d+1):
            search_sequences += distance(s,i)
    seq1 = [to_regex(s) for s in search_sequences]
    seq2 = [to_regex(rc(s)) for s in search_sequences]
    hit_number = 0
    # sort_results stores all hits before sorting/filtering.
    # NOT recommended for patterns with many hits
    if args.sort_results:
        hits_sense = list(set(flatten(
            [
                [
                    (i.span()[0],i.span()[0]+jl)
                    for i in re.finditer(lookahead(j),genome[chrom])
                ] for j,jl in zip(seq1,[
                    len(re.sub('\[.+?\]','1',l))
                    for l in seq1
                    ])
            ]
        )))
        hits_antisense = list(set(flatten(
            [
                [
                    (i.span()[0],i.span()[0]+jl)
                    for i in re.finditer(lookahead(j),genome[chrom])
                ] for j,jl in zip(seq2,[
                    len(re.sub('\[.+?\]','1',l))
                    for l in seq2
                    ])
            ]
        )))
        all_hits = sorted(
            list(
                zip(
                    hits_sense,
                    ['+']*len(hits_sense)
                    )
                ) + \
            list(
                zip(
                    hits_antisense,
                    ['-']*len(hits_antisense)
                    )
                )
            )
        for h,s in all_hits:
            line_to_write = '\t'.join(
                [
                    chrom,
                    str(h[0]),
                    str(h[1]),
                    chrom + '_' + str(hit_number),
                    '0',
                    s
                ]
            ) + '\n'
            if parallel:
                queue.put(line_to_write)
            else:
                queue.write(line_to_write)
            hit_number+=1
    # Iterative search performs filtering and output dynamically
    # Output will not be sorted
    else:
        start_end_pairs = {}
        for j,jl in zip(seq1,[len(re.sub('\[.+?\]','1',l)) for l in seq1]):
            pattern = re.compile(lookahead(j))
            for i in pattern.finditer(genome[chrom]):
                hit_span = (i.span()[0],i.span()[0]+jl)
                if not hit_span[0] in start_end_pairs:
                    start_end_pairs[hit_span[0]] = set([hit_span[1]])
                elif not hit_span[1] in start_end_pairs[hit_span[0]]:
                    start_end_pairs[hit_span[0]].add(hit_span[1])
                else:
                    continue
                hit_number += 1
                line_to_write = '\t'.join(
                    [
                        chrom,
                        str(hit_span[0]),
                        str(hit_span[1]),
                        chrom + '_' + str(hit_number),
                        '0',
                        '+'
                    ]
                ) + '\n'
                if parallel:
                    queue.put(line_to_write)
                else:
                    queue.write(line_to_write)
        
        start_end_pairs = {}
        for j,jl in zip(seq2,[len(re.sub('\[.+?\]','1',l)) for l in seq2]):
            pattern = re.compile(lookahead(j))
            for i in pattern.finditer(genome[chrom]):
                hit_span = (i.span()[0],i.span()[0]+jl)
                if not hit_span[0] in start_end_pairs:
                    start_end_pairs[hit_span[0]] = set([hit_span[1]])
                elif not hit_span[1] in start_end_pairs[hit_span[0]]:
                    start_end_pairs[hit_span[0]].add(hit_span[1])
                else:
                    continue
                hit_number += 1
                line_to_write = '\t'.join(
                    [
                        chrom,
                        str(hit_span[0]),
                        str(hit_span[1]),
                        chrom + '_' + str(hit_number),
                        '0',
                        '-'
                    ]
                ) + '\n'
                if parallel:
                    queue.put(line_to_write)
                else:
                    queue.write(line_to_write)
    
    print('{} hits: {}'.format(chrom,hit_number))
    

IUPAChash = {}
IUPACcomp = {}
IUPACnot  = {}
values = [
    'A','C','G','T','T','',
    '[AC]','[AG]','[AT]','[CG]','[CT]','[GT]',
    '[ACG]','[ACT]','[AGT]','[CGT]','.','.'
]
keys = [
    'A','C','G','T','U','-',
    'M','R','W','S','Y','K',
    'V','H','D','B','N','.'
]
complements = [
    'T','G','C','A','A','',
    'K','Y','W','S','R','M',
    'B','D','H','V','N','.'
]
notkeys = [
    'B','D','H','V','V','-',
    'K','Y','S','W','R','M',
    'T','G','C','A','N','.'
]

for k,v,c,n in zip(keys,values,complements,notkeys):
     IUPAChash[k] = v
     IUPACcomp[k] = c
     IUPACnot[k] = n

#################################
# LOADING DATA FROM INPUT FILES #
#################################

genome = fu.import_genome(args.genome_FASTA,keep_case=False)
chromosomes = {}

for k,v in genome.items():
    chromosomes[k] = len(v)

linecounter = 0
if __name__ == '__main__':
    for line in open(args.match_file):
        linecounter += 1
        if line[0] == '#':
            continue
        try:
            BED_filename,search_sequence,mismatches = line.rstrip().split(' ')
        except:
            print('Wrong number of columns in line '+str(linecounter))
            continue
        if args.parallel:
            queue = mp.Queue()
            STOP_TOKEN = "FINISHED"
            writer_process = mp.Process(
                target=writer,
                args=(out_folder+BED_filename+'.bed',queue,STOP_TOKEN)
            )
            writer_process.start()
            all_threads = []
            chromosomes_by_size = [
                k for v,k in sorted(
                    [(v,k) for k,v in chromosomes.items()]
                    ,reverse=True
                )
            ]
            print("SEARCH: {}".format(search_sequence))
            for chrom in chromosomes_by_size:
                all_threads.append(
                    mp.Process(
                        target=IUPAC_matches,
                        args=(
                            search_sequence,
                            mismatches,
                            chrom,
                            queue,
                            True
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
            print("SEARCH: {}".format(search_sequence))
            queue = open(out_folder+BED_filename+'.bed','w')
            for chrom in sorted(chromosomes.keys()):
                IUPAC_matches(
                    search_sequence,
                    mismatches,
                    chrom,
                    queue,
                    False
                )
            
            queue.close()

        
    
print('Sequence search complete.')
