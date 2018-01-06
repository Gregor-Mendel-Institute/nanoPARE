'''
Command-line: [1]merged_file [2]out1 [3]out2
'''
import sys
args = sys.argv[1:]
merged_file,out1,out2 = args

file          = open(merged_file)
split1        = open(out1,'w')
split2        = open(out2,'w')
entry         = []
linecounter   = 0
matepair      = 1
for line in file:
    line=line.rstrip()
    linecounter+=1
    if linecounter % 4 == 1:
        entry=[]
    entry.append(line)
    if linecounter % 4 == 0:
        if matepair == 1:
            split1.write('\n'.join(entry)+'\n')
            matepair = 2
        elif matepair == 2:
            split2.write('\n'.join(entry)+'\n')
            matepair = 1
split1.close()
split2.close()