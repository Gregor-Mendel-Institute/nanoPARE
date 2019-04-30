import sys
import os
import argparse
if sys.version_info >= (3,0):
    izip = zip
else:
    from itertools import izip

parser = argparse.ArgumentParser()
parser.add_argument(
    "--minlen", dest='MINLEN',
    help="Minimum length for a FASTQ read to keep.",
    default=12, type=int
)
parser.add_argument(
    "-O", "--output", dest='OUTPUT', type=str, 
    default=['cleaned.1.fastq','cleaned.2.fastq'], nargs=2,
    help="Filepaths to write processed FASTQ files."
)
parser.add_argument(
    "FILENAME", nargs=2
)
args = parser.parse_args()

file1=open(args.FILENAME[0])
file2=open(args.FILENAME[1])
tmp1=open(args.OUTPUT[0],'w')
tmp2=open(args.OUTPUT[1],'w')
print("Beginning to clean {} and {}...".format(args.FILENAME[0],args.FILENAME[1]))
entry1=[]
entry2=[]
linecounter=0
tooshortcount=0
for line1, line2 in izip(file1, file2):
    line1 = line1.rstrip()
    line2 = line2.rstrip()
    linecounter+=1
    if linecounter % 4 == 1:
        entry1=[]
        entry2=[]
    entry1.append(line1)
    entry2.append(line2)
    if linecounter % 4 == 0:
        if len(entry1[1]) >= args.MINLEN and len(entry2[1]) >= args.MINLEN:
            tmp1.write('\n'.join(entry1)+'\n')
            tmp2.write('\n'.join(entry2)+'\n')
        else:
            tooshortcount+=1
tmp1.close()
tmp2.close()
print("Finished cleaning.")
print(str(tooshortcount),"read pairs removed.")
