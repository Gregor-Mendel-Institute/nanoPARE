#sRNA_shuffler.py
#input:
##1. FASTA file of sRNAs for which sequences should be shuffled using the random() function
##2. Number of shuffled sets to be produced
#output
##1. X number of shuffled sRNA files

import sys, random

inFile = sys.argv[1].strip()
shuffNum = int(sys.argv[2].strip())

#read in .fa file
ifh = open(inFile)
inLines = ifh.readlines()
ifh.close()

for j in range(0,shuffNum):
    outFile = inFile.split('.fa')[0] + '.' + str(j) + '.shuffled.fa'
    ofh = open(outFile,'w')

    for i in range(0,len(inLines),2):
        entry = inLines[i].strip() + '_shuffled_' + str(j)
        seq = list(inLines[i + 1].strip())
        random.shuffle(seq)
        seq_shuffled = ''
        for base in seq:
            seq_shuffled = seq_shuffled + base
        print >>ofh,entry
        print >>ofh,seq_shuffled
    
ofh.close()
