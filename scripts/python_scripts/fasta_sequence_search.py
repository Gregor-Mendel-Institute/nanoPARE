'''
Commandline arguments: [1] genome_FASTA [2] match_file [3] out_folder
genome_FASTA - Full path of the FASTA file to search
match_file   - Full path of a space-delimited file with information on searches to perform.
    Each line must contain the following information:
        BED.filename: Path to a BED file that will be written containing all hits
        search.sequence(s): DNA sequence regex pattern that follows IUPAC nucleotide naming conventions.
            Separate multiple search sequences with a comma
        mismatches: Number of mismatches allowed in a matching sequence
            Separate multiple search sequences with a comma
'''
############################
# GLOBAL ENVIRONMENT SETUP #
############################

import sys
import re,os
genome_FASTA = sys.argv[1]
match_file   = sys.argv[2]
try:
    out_folder = sys.argv[3]
except:
    out_folder = './'

if out_folder[-1] != '/':
    out_folder+='/'

genome={}
chromosomes={}

def flatten(list_of_lists):
    return [item for sublist in list_of_lists for item in sublist]

def distance(sequencelist,dist,wildcard=False):
    if type(sequencelist) is str:
        sequencelist = [sequencelist]
    for i in range(dist):
        if wildcard:
            sequencelist = list(set(flatten([['.'.join([string[0:j],string[j+1:]]) for j in range(len(string))] for string in sequencelist])))
        else:
            sequencelist = list(set(flatten([[(IUPACnot[string[j]]).join([string[0:j],string[j+1:]]) for j in range(len(string))] for string in sequencelist])))
    return sequencelist

'''
Modifies a regex expression so that the first character is followed by a lookahead assertion for the remainder of the string
'''
def lookahead(string):
    if string[0] != '[':
        cutposition = 0
    else:
        cutposition = min([i for i in range(len(string)) if string[i]==']'])
    return string[0:cutposition+1]+'(?='+string[cutposition+1:]+')'

def rc(sequence):
    return ''.join(reversed([IUPACcomp[i] for i in sequence.upper()]))

def to_regex(sequence):
    return ''.join([IUPAChash[i] for i in sequence.upper()])

def genome_matches(sequence,dist=0,name='matches',quiet=False):
    sequence = sequence.split(',')
    seqlen   = [len(i) for i in sequence]
    if type(dist) is str:
        dist = [int(i) for i in dist.split(',')]
    elif type(dist) is int:
        dist = [dist]*len(sequence)
    search_sequences = []
    assert len(sequence) == len(dist)
    for s,d in zip(sequence,dist):
        print(s,d)
        search_sequences += distance(s,d)
    out_file=name+'.bed'
    out_bed=open(out_file,'w')
    seq1=[to_regex(s) for s in search_sequences]
    seq2=[to_regex(rc(s)) for s in search_sequences]
    hit_number=0
    for chrom in chromosomes.keys():
        hits_sense=list(set(flatten([[(i.span()[0],i.span()[0]+jl) for i in re.finditer(lookahead(j),genome[chrom])] for j,jl in zip(seq1,[len(re.sub('\[.+?\]','1',l)) for l in seq1])])))
        hits_antisense=list(set(flatten([[(i.span()[0],i.span()[0]+jl) for i in re.finditer(lookahead(j),genome[chrom])] for j,jl in zip(seq2,[len(re.sub('\[.+?\]','1',l)) for l in seq2])])))
        all_hits = sorted(list(zip(hits_sense,['+']*len(hits_sense)))+list(zip(hits_antisense,['-']*len(hits_antisense))))
        for h,s in all_hits:
            out_bed.write('\t'.join([chrom,str(h[0]),str(h[1]),'match'+str(hit_number),'0',s])+'\n')
            hit_number+=1
    out_bed.close()
    if not quiet:
        print('Total hits: '+str(hit_number))

IUPAChash = {}
IUPACcomp = {}
IUPACnot  = {}
values=['A','C','G','T','T','',
      '[AC]','[AG]','[AT]','[CG]','[CT]','[GT]',
      '[ACG]','[ACT]','[AGT]','[CGT]','.','.']
keys=['A','C','G','T','U','-',
        'M','R','W','S','Y','K',
        'V','H','D','B','N','.']
complements=['T','G','C','A','A','',
             'K','Y','W','S','R','M',
             'B','D','H','V','N','.']
notkeys=['B','D','H','V','V','-',
        'K','Y','S','W','R','M',
        'T','G','C','A','N','.']
for k,v,c,n in zip(keys,values,complements,notkeys):
     IUPAChash[k]=v
     IUPACcomp[k]=c
     IUPACnot[k]=n

#################################
# LOADING DATA FROM INPUT FILES #
#################################

chrom='none'
genome_file=open(genome_FASTA)
for line in genome_file:
    line=line.rstrip()
    if line[0]=='>':
        if chrom!='none':
            genome[chrom]=''.join(x)
        chrom=line[1:len(line)]
        x=[]
        continue
    x.append(line)
genome[chrom]=''.join(x)

for k,v in genome.items():
    chromosomes[k]=len(v)

linecounter=0
for line in open(match_file):
    linecounter+=1
    if line[0]=='#':
        continue
    try:
        BED_filename,search_sequence,mismatches = line.rstrip().split(' ')
    except:
        print('Wrong number of columns in line '+str(linecounter))
        continue
    genome_matches(search_sequence,mismatches,out_folder+BED_filename)
    


