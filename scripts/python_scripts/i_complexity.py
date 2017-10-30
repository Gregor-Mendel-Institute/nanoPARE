# -*- coding: utf-8 -*-

# As defined in Becher & Heiber (2012) Theoretical Computer Science: "A linearly computable measure of string complexity"
# Forward Repetitions Array = F
# Backward Repititions Array = B
# x^k defined by function K
# N = number of unique substrings of length l
# M = number of elements in Bs < l
# doms tells if x dominates y

# SRA = Sorted Repititions Array
# Z = Sorted Repetitions Array of Extended de Bruijn String
# C = Most dominant Sorted Repetitions Array of length l
##############
# For every sequence s in A*, SRA(s) is bounded by Z(len(s),a) and C(len(s))
# [doms(SRA(s),Z(len(s),a)) and doms(C(len(s)),SRA(s)) for all s in A*]
##############

# S = Suffix array of sequence s
# LCP = Longest Common Prefix of sequence s
# prefixes = an array of sorted LCP prefixes
# dlog = Discrete derivative function of log

# Ia = I-complexity function for any repetitions array A of alphabet a
# Is = I-complexity function for any string s of alphabet a

import math,os
from collections import defaultdict

def K(x,k):
    sl=k*len(x)
    return x*int(sl//len(x))+x[0:int(sl%len(x))]

def doms(x,y):
    return all([a>=b for a,b in zip(x,y)])
    
def rev(x):
    return x[::-1]

def F(s):
    Fs=[0]*len(s)
    for i in range(len(s)):
        Lsuffixes=[s[i:j+1] for j in range(i,len(s))]
        Ssub=s[(i+1):len(s)]
        Fs[i]=max([len(x) for x in Lsuffixes if x in Ssub]+[0])

    return Fs

def B(s):
    Bs=[0]*len(s)
    for i in range(len(s)+1):
        Lprefixes=[s[i-(j+1):i] for j in list(range(i))]
        Ssub=s[0:(i-1)]
        Bs[i-1]=max([len(x) for x in Lprefixes if x in Ssub]+[0])
    return Bs

def N(s,l):
    return len(set([s[i:(i+l)] for i in range(len(s)-l+1)]))

def M(array,l):
    return len([j for j in [i<l for i in array] if j])

def de_bruijn(a, k, wrap=True):      # De Bruijn sequence for alphabet a and subsequences of length k.
    try:
        _ = int(a)
        alphabet = list(map(str, range(a)))
    except (ValueError, TypeError):
        alphabet = a
        a = len(a)
    b = [0] * a * k
    sequence = []
    def db(t, p):
        if t > k:
            if k % p == 0:
                sequence.extend(b[1:p + 1])
        else:
            b[t] = b[t - p]
            db(t + 1, p)
            for j in range(b[t - p] + 1, a):
                b[t] = j
                db(t + 1, t)
    db(1, 1)
    dbs = "".join(alphabet[i] for i in sequence)
    if wrap: dbs = dbs + dbs[0:k-1]
    return dbs

def SRA(s):
    Bs=B(s)
    Bs.sort()
    return Bs

def Z(l,a):
    Zl=[0]*l
    k=0
    for i in range(l):
        if i>=a**(k+1)+k:
            k+=1
        Zl[i]=k
    return Zl

def C(l):
    return list(range(l))

def S(s):
    j=[0]*len(s)
    for i in range(len(s)):
        j[i]=s[i:len(s)]
    k=sorted(j)
    Ss=[0]*len(s)
    for i in range(len(k)):
        Ss[i]=[a for a,b in enumerate(j) if b==k[i]][0]
    return Ss

def sort_bucket(str, bucket, order): 
    d = defaultdict(list) 
    for i in bucket: 
        key = str[i:i+order] 
        d[key].append(i) 
    result = [] 
    for k,v in sorted(d.items()): 
        if len(v) > 1: 
            result += sort_bucket(str, v, order*2) 
        else: 
            result.append(v[0]) 
    return result 

def SMM(str): #suffix_array_ManberMyers from Ben Fulton @ http://algorithmicalley.com/archive/2013/06/30/suffix-arrays.aspx
    return sort_bucket(str, [i for i in range(len(str))], 1) 

def LCP(s):
    lcp = [0]*(len(s)-1)
    sa = SMM(s)
    for i in range(len(lcp)):
        l=1
        while s[sa[i]:sa[i]+l] == s[sa[i+1]:sa[i+1]+l]: l+=1
        lcp[i]=l-1
    return lcp

def LCPsa(s,sa):
    lcp = [0]*(len(s)-1)
    for i in range(len(lcp)):
        l=1
        while s[sa[i]:sa[i]+l] == s[sa[i+1]:sa[i+1]+l]: l+=1
        lcp[i]=l-1
    return lcp

assert LCPsa('mississippi',SMM("mississippi")) == LCP("mississippi")

# iSA and iLCP are based on the induced sorting algorithms in Fischer 2011
def SL(s): #Demarcates the S-suffix type at each position. If type S, then 1; elif type S*, then 2; else 0.
    n=len(s)
    sl=[0]*n
    sl[n-1]=1
    for i in reversed(range(n-1)):
        if s[i]<s[i+1] or (s[i]==s[i+1] and sl[i+1]): sl[i]=1
    for i in range(1,n):
        if sl[i] and sl[i-1]==0: sl[i]=2
    return sl

def S_(s,sl): #Produces an array of S*-suffixes
    R=[]
    start=False
    stop=False
    for i in range(len(s)):
        if sl[i]==2:
            if not start: start=i
            elif not stop:
                stop=i
                R.append(s[start:stop+1])
                start=stop
                stop=False
    return R

def iSA(s): #induced Suffix Array; output can be used as a suffix array
    sl=SL(s)
    R=S_(s,sl)
    names=[0]*len(R)
    SR=[(k,v) for v,k in enumerate(R)]
    order=dict([(v,k+1) for k,v in enumerate(sorted(set(R)))])
    for v,i in SR:
        names[i]=order[v]
    print(names)
    if len(set(R))<len(R):
        print(chash.to_hash(names))
        names=iSA(chash.to_hash(names))
    SA_=names
    return SA_

    
########################
def dlog(x,n=math.e):
    return math.log(x+1,n)-math.log(x,n)

def minalpha(s):
    return len(set(list(s)))

def Ia(A,a):
    return round(sum([dlog(i+1,a) for i in A]),4)

def Is(s,a='default'):
    if a=='default':a=minalpha(s)
    return Ia([0]+LCP(s),a)

##
def suffixes(s,SA):
    return 0

def prefixes(s):
    Ss=SMM(s)
    lcp=[0]+LCPsa(s,Ss)
    prf=[0]*len(lcp)
    for i in range(len(lcp)):
        prf[i]=s[Ss[i]:Ss[i]+lcp[i]]
    return prf

def allprefixes(s):
    Ss=SMM(s)
    lcp=[0]+LCPsa(s,Ss)
    prf=[0]*len(lcp)
    for i in range(len(lcp)):
        prf[i]=s[Ss[i]:Ss[i]+lcp[i]]
    return sorted(list(prf),key=len,reverse=False)

def LUSTR(s):
    prf=prefixes(s)
    frp=sorted([i[::-1] for i in prf])
    j=[frp[i] not in frp[i+1] for i in range(len(frp)-1)]+[True]
    prf= sorted(list(set([a[::-1] for a,b in zip(frp,j) if b])),key=len,reverse=True)
    return prf
    j=[1]*len(prf)
    for i in range(len(prf)):
        count=0
        for cur in prf[0:i+1]:
            if prf[i] in cur:
                count+=1
                if count>1:
                    j[i]=0
                    continue
    y=sorted(list(set([a for a,b in zip(prf,j) if b])),key=len,reverse=True)
    return y

def LUSTR_prebuilt(s,SA,LCP):
    prf=[0]*len(LCP)
    for i in range(len(LCP)):
        prf[i]=s[SA[i]:SA[i]+LCP[i]]
    prf=sorted(prf,key=len,reverse=True)

def mode(s,which="F"):
    if which=="F":
        Fs=F(s)
        return [s[l:l+max(Fs)] for l in [i for i,j in enumerate(Fs) if j==max(Fs)]]
    elif which=="B":
        Fs=B(s)
        return [s[l-max(Fs)+1:l+1] for l in [i for i,j in enumerate(Fs) if j==max(Fs)]]

def simplify(s,keepnums=True):
    if keepnums:
        dic1=".':,-;_!/\"\n\t*()?{}[]#@%$"
    else:
        dic1=".':,-;_!/\"\n\t*()?{}[]#@%$0123456789"
    dic2=" "*len(dic1)
    return s.translate(s.maketrans(dic1,dic2)).replace(" ","").lower()

def ACGTN(string):
    dic1="WSKDMRY"
    dic2="N"*len(dic1)
    return string.translate(string.maketrans(dic1,dic2)).replace(" ","")
##

 
s="A THIEF IS A THING THAT A THISTLE FARMER HATES, BUT A THIEF IS ALSO A TIRESOME FARMHAND"

Fs=F(s)
Fsr=F(rev(s))
Bs=B(s)
Bsr=B(rev(s))
Fs.sort()
Fsr.sort()
Bs.sort()
Bsr.sort()

assert Fs==Fsr and Fs==Bs and Fs==Bsr #THEOREM 9
assert doms(B("FARMER FARNSWORTH, "+"FRANNY FORTUNE"),B("FARMER FARNSWORTH, ")+B("FRANNY FORTUNE")) #LEMMA 11a
assert not doms(B("FARMER FARNSWORTH, ")+B("FRANNY FORTUNE"),B("FARMER FARNSWORTH, "+"FRANNY FORTUNE")) #LEMMA 11a
assert B("FARMER FARNSWORTH,")+B("BILL_BLIZZY") == B("FARMER FARNSWORTH,"+"BILL_BLIZZY")    #LEMMA 11b
assert [N(s,l)+l-2 for l in range(1,13)] == [M(LCP(s),l) for l in range(1,13)]  #LEMMA 17
assert SRA(s) == [0]+sorted(LCP(s)) #THEOREM 18
assert max([Is("FARMERFARNSWORTH",26),Is("FRANNYFORTUNE",26)]) <= Is("FARMERFARNSWORTHFRANNYFORTUNE",26)    #THEOREM 25.3a
assert Is("FARMERFARNSWORTHFRANNYFORTUNE",26) <= Is("FARMERFARNSWORTH",26)+Is("FRANNYFORTUNE",26)   #THEOREM 25.3b
assert Is("BIBBLEBOBBLE",26) <= Is("EBIBBLEBOBBLE",26) and Is("EBIBBLEBOBBLE",26) <= Is("BIBBLEBOBBLE",26)+1 #THEOREM 25.4a
assert Is("BIBBLEBOBBLE",26) <= Is("BIBBLEBOBBLEB",26) and Is("BIBBLEBOBBLEB",26) <= Is("BIBBLEBOBBLE",26)+1 #THEOREM 25.4b
assert Ia(SRA(s),28)==Is(s,28)
assert Ia(B(s),28)==Ia(sorted(B(s)),28)

def windowsplit(string,order):
    return [string[i-order:i+order+1] for i in range(order,len(string)-order)]


