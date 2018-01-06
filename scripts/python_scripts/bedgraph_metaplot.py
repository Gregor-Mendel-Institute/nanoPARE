import re,sys

usage="\
Outputs a column of count values representing the aggregate value in a BEDGRAPH file of all features in a GFF file\
\n\
Commandline arguments:\n\
    -P=[input plus bedgraph]      (required filepath)\n\
    -M=[input minus bedgraph]      (required filepath)\n\
    -G=[input GFF]                (required filepath)\n\
    -L=[lengths table]            (required filepath; output from fasta_lengths.py)\n\
    -S=[strand]                   (default: both; options: both|. , plus|+ , minus|-)\n\
    -F=[feature type]             (default: exon; options: any string)\n\
    -V=[vertical normalization]   (default: norm; options: max,mean,norm)\n\
    -H=[horizontal normalization] (default: left; options: left,right,center,stretch<int>)\n\
    -M=[max length]               (default: 1000; options: int)\n\
    -A=[add length]               (default: 0; options: int)\n\
    -splice=[splice out introns]  (default: True; options: bool)\n\
\n\
This vector of values will grow in size over time based on the rule provided in -H:\n\
    left:         left-aligned\n\
    right:        right-aligned\n\
    center:       center-aligned\n\
    stretch<int>: All features will be resized into a vector of fixed length <int>\n\
\n\
Values added to the vector will be normalized prior according to -V:\n\
    max:          Sum of all values\n\
    mean:         Mean of all values\n\
    norm:         The values in each feature are scaled / max(abs(values))\n\
"

STRAND     = 'both'
FEATURE    = 'exon'
VERTICAL   = 'norm'
HORIZONTAL = 'left'
MAXLEN     = 1000
ADDLEN     = 0
SPLICE     = True

if len(sys.argv) < 2:
    print(usage)
    sys.exit()
args = sys.argv[1:]
for arg in args:
    if len(arg.split('=')) == 1:
        print("Please connect arguments and values with '='")
        print(usage)
        sys.exit()
    option,value=arg.split('=')
    if option == '-P':   PLUS_BG = value
    if option == '-M':   MINUS_BG = value
    elif option == '-G': GFF_IN      = value
    elif option == '-L': LENGTHS     = value
    elif option == '-S': STRAND      = value
    elif option == '-F': FEATURE     = value
    elif option == '-V': VERTICAL    = value
    elif option == '-H': HORIZONTAL  = value
    elif option == '-splice': SPLICE = bool(value)
    elif option == '-maxlen': MAXLEN      = int(value)
    elif option == '-A': ADDLEN      = int(value)
if 'PLUS_BG' not in globals():
    print("ERROR: missing required argument -P=[input bedgraph]")
    print(usage)
    sys.exit()
if 'MINUS_BG' not in globals():
    print("ERROR: missing required argument -M=[input bedgraph]")
    print(usage)
    sys.exit()
if 'GFF_IN' not in globals():
    print("ERROR: missing required argument -G=[input GFF]")
    print(usage)
    sys.exit()
if 'LENGTHS' not in globals():
    print("ERROR: missing required argument -L=[lengths table]")
    print(usage)
    sys.exit()

def which(x,value=True):
    """Returns a list of locations in x that satisfy value"""
    return [a for a,b in enumerate(x) if b==value]


def notwhich(x,value=0):
    """Returns a list of locations in x that do not satisty value"""
    return [a for a,b in enumerate(x) if b!=value]


def flatten(list_of_lists):
    """Collapses a list/tuple of lists into a single list"""
    return [item for sublist in list_of_lists for item in sublist]
    
'''
'chromosomes' contains the lengths of all chromosomes the that BEDGRAPH contains values for.
Expects a two-column tab-separated file with:
    chromosome  length
Provided with the -L argument.
''' 
chromosomes={}
lengths_file=open(LENGTHS)
for line in lengths_file:
    chrom,length=line.rstrip().split('\t')
    chromosomes[chrom]=int(length)

'''
'coverage' is a dictionary of float vectors for each nucleotide in the genome.
Contains the value of the BEDGRAPH file at each position.
BEDGRAPH format is zero-based, half-open
'''
coverage_plus={}
for chrom,chromlen in chromosomes.items():
    coverage_plus[chrom] = [0]*chromlen
coverage_file = open(PLUS_BG)

for line in coverage_file:
    chrom,start,end,count = line.rstrip().split()
    count = float(count)
    for i in range(int(start),int(end)):
        coverage_plus[chrom][i] = count

coverage_minus={}
for chrom,chromlen in chromosomes.items():
    coverage_minus[chrom] = [0]*chromlen

coverage_file = open(MINUS_BG)
for line in coverage_file:
    chrom,start,end,count = line.rstrip().split()
    count = float(count)
    for i in range(int(start),int(end)):
        coverage_minus[chrom][i] = count


if VERTICAL.lower() == 'rpm':
    million_counts = sum(
        [
            sum(coverage_plus[i])
            for i in chromosomes.keys()
        ] + [
            sum(coverage_minus[i])
            for i in chromosomes.keys()
        ]
    )*(10**-6)

'''
Begins reading the input GFF file and extracting the features of type -F.
'value_vector' is the set of running values that are aggregated.
The vector will grow in size over time based on the rule provided in -H:
    left:         left-aligned
    right:        right-aligned
    center:       center-aligned
    stretch<int>: All features will be resized into a vector of fixed length <int>
    
Values added to the vector will be normalized prior according to -V:
    max:          Sum of all values
    mean:         Mean of all values
    norm:         The values in each feature are scaled / max(abs(values))
'''
value_vector=[]
featurecount=0
if 'stretch' in HORIZONTAL:
    HORIZONTAL = int(HORIZONTAL.replace('stretch',''))
    value_vector = [float(0)]*HORIZONTAL
    
# TODO: retool gff_utils to import gff3 format files
ref_transcripts = {}
refGFF = open(GFF_IN)
for line in refGFF:
    if line[0] == '#':continue
    chrom,source,gtype,start,end,score,strand,phase,other = line.rstrip().split('\t')
    if gtype != FEATURE:
        continue
    
    parent = re.search('^.*;?Parent=([^;,]+),?;?.*$',other)
    if parent:
        parent = parent.groups()[0]
    ID = re.search('(transcript_id|ID)=(.+?)[\.:;].+$',other)
    if ID:
        ID = ID.groups()[1]
    else:
        ID = parent.split('\.')[0]
    
    ID = ID.replace('-Protein','')
    iso_length = len(re.findall(',',parent))+1
    isos = re.search(','.join([ID+'(\.[0-9]+)?']*iso_length),parent).groups()
    iso_ID = parent.split(',')
    iso_ID = [i for i in iso_ID if 'Protein' not in i]
    for i in iso_ID:
        ref_transcripts[i] = ref_transcripts.get(i,{})
        if gtype != 'miRNA_primary_transcript':
            exon_number = re.search(
                'ID=' + ID + ':(pseudogenic_)?exon:([0-9]+).*;Parent.+$',
                other
            )
            if exon_number:
                exon_number = exon_number.groups()[1]
            else:
                exon_number = str(max([int(j) for j in ref_transcripts[i].get('exon_nums',['0'])])+1)
        else:
            exon_number = '1'
    
    for i in isos:
        if i:
            iso_ID = ID + i
        else:
            iso_ID = ID
        
        ref_transcripts[iso_ID]['chromosome'] = chrom
        ref_transcripts[iso_ID]['strand'] = strand
        ref_transcripts[iso_ID]['start'] = ref_transcripts[iso_ID].get('start',[])+[int(start)]
        ref_transcripts[iso_ID]['end'] = ref_transcripts[iso_ID].get('end',[])+[int(end)]
        ref_transcripts[iso_ID]['exon_nums'] = ref_transcripts[iso_ID].get('exon_nums',[])+[exon_number]  
        
ref_IDs = sorted(list(ref_transcripts.keys()))
print('# {} reference transcripts: {}'.format(
    len(ref_IDs),
    GFF_IN
))


for ID in ref_IDs:
    transcript = ref_transcripts[ID]
    
    chrom = transcript['chromosome']
    # load 'start' and 'end' files listed in the GFF,
    # adjusting for 1-based closed indexing by subtracting 1
    # from the start location
    start = min(transcript['start']) - 1 - ADDLEN
    end = max(transcript['end']) + ADDLEN
    strand = transcript['strand']
    if chrom not in chromosomes:
        continue
    
    if start < 0 or end > chromosomes[chrom]:
        continue
    
    if strand == '+' and STRAND.lower() in ['-','minus','m']:
        continue
    
    if strand == '-' and STRAND.lower() in ['+','plus','p']:
        continue
    
    if SPLICE:
        exon_starts = sorted([int(i) for i in transcript['start']])
        exon_ends = sorted([int(i) for i in transcript['end']])
        exon_starts[0] = exon_starts[0] - ADDLEN
        exon_ends[-1] = exon_ends[-1] + ADDLEN
        positions = flatten(
            [list(range(a - 1, b)) for a,b in zip(exon_starts, exon_ends)]
        )
        if strand == '+':
            line_values = [coverage_plus[chrom][i] for i in positions]
        elif strand == '-':
            line_values = [coverage_minus[chrom][i] for i in positions]
        # print(len(positions))
    else:
        if strand == '+':
            line_values = coverage_plus[chrom][start:end]
        elif strand == '-':
            line_values = coverage_plus[chrom][start:end]
    
    if strand == '-':
        line_values = list(reversed(line_values))
    
    if VERTICAL == 'norm':
        norm_value  = float(max([abs(i) for i in line_values]))
        if norm_value == 0:
            continue
        line_values = [float(i)/norm_value for i in line_values]
    
    if type(HORIZONTAL) is int:
        bins = [float(0)]*HORIZONTAL
        partition = float(len(line_values))/HORIZONTAL
        for i in range(0,HORIZONTAL):
            values_to_bin = line_values[int(round(partition*(i))):int(round(partition*(i+1)))]
            if not len(values_to_bin)==0:
                if VERTICAL == 'max':
                    bins[i] = sum(values_to_bin)
                else:
                    bins[i] = sum(values_to_bin)/len(values_to_bin)
        value_vector = [a+b for a,b in zip(value_vector,bins)]
    else:
        veclen = len(value_vector)
        newlen = len(line_values)
        if HORIZONTAL == 'left':
            if veclen >= newlen:
                value_vector = [a+b for a,b in zip(value_vector,line_values+[0]*(veclen-newlen))][:MAXLEN]
            else:
                value_vector = [a+b for a,b in zip(value_vector+[0]*(newlen-veclen),line_values)][:MAXLEN]
        elif HORIZONTAL == 'right':
            if veclen >= newlen:
                value_vector = [a+b for a,b in zip(value_vector,[0]*(veclen-newlen)+line_values)][-MAXLEN:]
            else:
                value_vector = [a+b for a,b in zip([0]*(newlen-veclen)+value_vector,line_values)][-MAXLEN:]
        elif HORIZONTAL == 'center':
            print('Error: center-align not yet implemented')
            sys.exit()
    featurecount+=1

if VERTICAL in ['mean','norm']:
    value_vector = [i/featurecount for i in value_vector]

if VERTICAL.lower() == 'rpm':
    value_vector = [i/million_counts for i in value_vector]

for v in value_vector:
    print('{:.10f}'.format(v))








