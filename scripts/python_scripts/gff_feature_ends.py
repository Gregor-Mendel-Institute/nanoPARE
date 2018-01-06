
'''
Commandline arguments: -I=[input file] -O=[output file] -E=[end type] -D=[distance] -F=[feature type] -T=[training transcripts]
Defaults:
    -I NONE (required)
    -O NONE (required)
    -E TSS  (options: TSS|start|s|S|5|5' ; TES|TTS|e|E|end|stop|3|3')
    -D 1000 (options: int)
    -F gene (options: string)
'''
import re,sys

ENDTYPE  = 'TSS'
DISTANCE = 1000
FEATURE  = 'mRNA'
TRAINING = 'none'

args = sys.argv[1:]
for arg in args:
    option,value=arg.split('=')
    if option == '-I':
        GFF_IN=value
    elif option == '-O':
        GFF_OUT=value
    elif option == '-E':
        if value in ['TSS','start','s','S','5',"5'"]:
            ENDTYPE='TSS'
        elif value in ['TES','TTS','e','E','end','stop','3',"3'"]:
            ENDTYPE='TES'
    elif option == '-D':
        DISTANCE=int(value)
    elif option == '-F':
        FEATURE=value
    elif option == '-T':
        TRAINING=value
    
GFF=open(GFF_IN)
GFFout=open(GFF_OUT,'w')
if TRAINING != 'none':
    use_all_genes=False
    training_data=set([l.rstrip() for l in open(TRAINING)])
else:
    use_all_genes=True
for line in GFF:
    if line[0] == '#':
        continue
    l=line.rstrip().split('\t')
    if l[2] != FEATURE:
        continue
    if not use_all_genes:
        try:
            ID = re.search('ID=([^;$]+);?.+$',l[8]).groups()[0]
        except:
            continue
        if ID not in training_data:
            continue
    start  = int(l[3])
    end    = int(l[4])
    strand = l[6]
    if strand == '+':
        if ENDTYPE == 'TSS':
            newcoords = [start-DISTANCE,start+DISTANCE]
        elif ENDTYPE == 'TES':
            newcoords = [end-DISTANCE,end+DISTANCE]
    elif strand == '-':
        if ENDTYPE == 'TSS':
            newcoords = [end-DISTANCE,end+DISTANCE]
        elif ENDTYPE == 'TES':
            newcoords = [start-DISTANCE,start+DISTANCE]
    else:
        continue
    if any([i<1 for i in newcoords]):
        continue
    newline='\t'.join(l[0:3]+[str(i) for i in newcoords]+l[5:])+'\n'
    GFFout.write(newline)
GFFout.close()
