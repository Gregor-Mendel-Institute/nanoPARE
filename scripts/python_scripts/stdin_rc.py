'''
Input arguments: [1]seq_in [2]complement_type ('RC','R','C'; default="RC")
'''
import sys
seq_in          = sys.argv[1]

try:
    complement_type = sys.argv[2]
except:
    complement_type = 'RC'
    
complement_type = complement_type.upper()
patternhash={'A':'T','T':'A','G':'C','C':'G','U':'A','R':'Y','Y':'R','S':'W','W':'S','K':'M','M':'K','B':'V','D':'H','H':'D','V':'B','N':'N','[':']',']':'[','-':'-','X':'X','{':'}','}':'{'}

def comp(seq):
    seq=seq.upper()
    return ''.join([patternhash[i] for i in list(seq)])

def rev(string):
    return ''.join(list(reversed(string)))
    
def rc(seq):
    seq=seq.upper()
    return ''.join(list(reversed([patternhash[i] for i in list(seq)])))

if complement_type == 'C':
    print(comp(seq_in))
elif complement_type == 'R':
    print(rev(seq_in))
elif complement_type in ['RC','CR']:
    print(rc(seq_in))