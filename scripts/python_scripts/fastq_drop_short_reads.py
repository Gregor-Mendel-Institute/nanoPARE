'''
Command-line: [1]directory [2]mate1 [3]out1 [4]minlen
'''
import sys
directory,mate1,out1,minlen=sys.argv[1:5]
minlen=int(minlen)
### USER DEFINED ###
# minlen=5
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
for i in file1:
    linecounter+=1
    line1=i.rstrip()
    if linecounter % 4 == 1:
        entry1=[]
    entry1.append(line1)
    if linecounter % 4 == 0:
        if len(entry1[1])>=minlen:
            tmp1.write('\n'.join(entry1)+'\n')
        else:
            tooshortcount+=1
tmp1.close()
print "Finished cleaning."
print str(tooshortcount),"reads removed."