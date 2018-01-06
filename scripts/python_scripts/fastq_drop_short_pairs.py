'''
Command-line: [1]directory [2]mate1 [3]mate2 [4]out1 [5]out2 [6]minlen
'''
import sys
directory,mate1,mate2,out1,out2,minlen=sys.argv[1:7]
minlen=int(minlen)
### USER DEFINED ###
# minlen=5
# directory='C:/Users/schon.admin/Desktop/'
# mate1='testfile1.txt'
# mate2='testfile2.txt'
# out1='testfile1_drop.txt'
# out2='testfile2_drop.txt'
####################
if not directory.endswith('/'):
    directory=directory+'/'
file1=open(directory+mate1)
file2=open(directory+mate2)
num_lines1 = sum(1 for line in file1)
num_lines2 = sum(1 for line in file2)
assert num_lines1 == num_lines2
print num_lines1,num_lines2
file1.close()
file2.close()
file1=open(directory+mate1)
file2=open(directory+mate2)
tmp1=open(directory+out1,'w')
tmp2=open(directory+out2,'w')
print "Beginning to clean..."
entry1=[]
entry2=[]
linecounter=0
tooshortcount=0
for i in range(num_lines1):
    linecounter+=1
    line1=file1.readline().rstrip()
    line2=file2.readline().rstrip()
    if linecounter % 4 == 1:
        entry1=[]
        entry2=[]
    entry1.append(line1)
    entry2.append(line2)
    if linecounter % 4 == 0:
        if len(entry1[1])>=minlen and len(entry2[1])>=minlen:
            tmp1.write('\n'.join(entry1)+'\n')
            tmp2.write('\n'.join(entry2)+'\n')
        else:
            tooshortcount+=1
tmp1.close()
tmp2.close()
print "Finished cleaning."
print str(tooshortcount),"read pairs removed."