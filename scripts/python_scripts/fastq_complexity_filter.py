# -*- coding: utf-8 -*-

'''
Command-line: [1]directory [2]mate1 [3]out1 [4]mincomp
'''
from i_complexity import *
import sys
directory,mate1,out1,mincomp=sys.argv[1:5]
mincomp=float(mincomp)
### USER DEFINED ###
# mincomp=5
# directory='C:/Users/schon.admin/Desktop/'
# mate1='testfile1.txt'
# out1='testfile1_drop.txt'
####################
if not directory.endswith('/'):
    directory=directory+'/'
# file1=open(directory+mate1)
# num_lines1 = sum(1 for line in file1)
# print str(num_lines1/4)
# print "Beginning to clean..."
# file1.close()
file1=open(directory+mate1)
tmp1=open(directory+out1,'w')
entry1=[]
linecounter=0
tooshortcount=0
table_out = True
if table_out:
    comptable = {}
for i in file1:
    linecounter+=1
    # line1=file1.readline().rstrip()
    line1=i.rstrip()
    if linecounter % 4 == 1:
        entry1=[]
    entry1.append(line1)
    if linecounter % 4 == 0:
        if len(entry1[1])>0:
            comp = Is(entry1[1],5)/len(entry1[1])
            if table_out:
                rcomp = round(comp,4)
                comptable[rcomp] = comptable.get(rcomp,0)+1
            if comp>=mincomp:
                tmp1.write('\n'.join(entry1)+'\n')
            else:
                tooshortcount+=1
        else:
            tooshortcount+=1
tmp1.close()
print "\ni-complexity filter: {}.".format(directory+mate1)
print str(tooshortcount),"low complexity reads removed.\n"
if table_out:
    outfile = open('comptable.tsv','w')
    for k,v in sorted([(k,v) for k,v in comptable.items()]):
        outfile.write('{}\t{}\n'.format(k,v))
    outfile.close()
