import os,re,sys

rc = {}
nucs = ['A','C','G','T','N']
cmps = ['T','G','C','A','N']

for n,c in zip(nucs,cmps):
    rc[n] = c
    
infile = open(sys.argv[1])
outfile = open(sys.argv[2],'w')
try:
    untemp_nuc = sys.argv[3]
except:
    untemp_nuc = "G"

untemp = {}
total = 0

for line in infile:
    if line[0]=='@':
        outfile.write(line)
        continue
    clipped_nuc = 'X'
    total += 1
    l = line.rstrip().split()
    SAMflags = bin(int(l[1]))[2:]
    SAMflags = '0'*(12-len(SAMflags))+SAMflags
    if SAMflags[-1] == '0' or SAMflags[-7] == '1':
        if SAMflags[-5] == '0':
            try:
                softclip = int(re.findall('(^[0-9]+)S',l[5])[0])
            except:
                softclip = 0
            if softclip != 0:
                clipped_nuc = l[9][(softclip - 1) : (softclip + 1)]
                untemp[clipped_nuc] = untemp.get(clipped_nuc,0) + 1
        elif SAMflags[-5] == '1':
            try:
                softclip = int(re.findall('([0-9]+)S$',l[5])[0])
            except:
                softclip = 0
            if softclip != 0:
                if softclip == 1:
                    clipped_nuc = ''.join(reversed([rc[i] for i in l[9][(-softclip -1) : ]]))
                else:
                    clipped_nuc = ''.join(reversed([rc[i] for i in l[9][(-softclip -1) : (-softclip +1)]]))
                untemp[clipped_nuc] = untemp.get(clipped_nuc,0) + 1
    if clipped_nuc[0] == 'A' and untemp_nuc == 'A': outfile.write(line)
    if clipped_nuc[0] == 'C' and untemp_nuc == 'C': outfile.write(line)
    if clipped_nuc[0] == 'G' and untemp_nuc == 'G': outfile.write(line)
    if clipped_nuc[0] == 'T' and untemp_nuc == 'T': outfile.write(line)
    if clipped_nuc[0] == 'X' and untemp_nuc == 'X': outfile.write(line)
outfile.close()

print(total)
patterns = sorted([(v,k) for k,v in untemp.items()],reverse=True)
for v,k in patterns:
    print(k+'\t'+str(v))
