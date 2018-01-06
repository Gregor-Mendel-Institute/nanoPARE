import sys
'''
Outputs a two-column tab-delimited file from a FASTA file:
    Feature name    Length of the feature sequence

Commandline arguments: [1] genome_FASTA
genome_FASTA - Full path of the FASTA file to search
'''

genome_FASTA = sys.argv[1]
split_on=' '

chrom='none'
genome={}
chromosomes={}
genome_file=open(genome_FASTA)
for line in genome_file:
    line=line.rstrip()
    if line[0]=='>':
        if chrom!='none':
            genome[chrom]=''.join(x)
        chrom=line[1:len(line)].split(split_on)[0]
        x=[]
        continue
    x.append(line)
if chrom!='none':
    genome[chrom]=''.join(x)

for k,v in genome.items():
    chromosomes[k]=len(v)

for k,v in sorted(list(chromosomes.items())):
    print('\t'.join([k,str(v)]))
