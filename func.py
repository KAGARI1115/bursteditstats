import numpy as np
import math
from math import comb,log,factorial,exp,ceil
from os import listdir
from os.path import join
from functools import lru_cache
import sys
from os.path import isfile, join
sys.setrecursionlimit(1000000)

def align(u,v):
    dp=[[0]*(len(v)+1) for i in range(len(u)+1)]
    arrow=[[0]*(len(v)+1) for i in range(len(u)+1)]
    for j in range(1,len(v)+1):
        dp[0][j]=dp[0][j-1]-1
        arrow[0][j]='right'
    for i in range(1,len(u)+1):
        dp[i][0]=dp[i-1][0]-1
        arrow[i][0]='down'
    for i in range(1,len(u)+1):
        for j in range(1,len(v)+1):
            downscore=dp[i-1][j]-1
            rightscore=dp[i][j-1]-1
            diagscore=dp[i-1][j-1]-int(u[i-1]!=v[j-1])
            if diagscore>=downscore and diagscore>=rightscore:
                arrow[i][j]='diag'
                dp[i][j]=diagscore
            elif downscore>=rightscore:
                arrow[i][j]='down'
                dp[i][j]=downscore
            else:
                arrow[i][j]='right'
                dp[i][j]=rightscore
    i,j=len(u),len(v)
    alignment=""
    while i>0 or j>0:
        direct=arrow[i][j]
        if direct=='diag':
            if u[i-1]==v[j-1]:
                alignment="m"+alignment
            else:
                alignment="*"+alignment
            i-=1
            j-=1
        elif direct=='down':
            alignment='/'+alignment
            i-=1
        else:
            alignment='-'+alignment
            j-=1
    return alignment,dp[-1][-1]

def listfiles(dir):
    ans = []
    for f in listdir(dir):
        newpath = join(dir,f) 
        if isfile(newpath):
            ans.append(newpath[10:])
        else:
            ans.extend(listfiles(newpath))
    return ans

def filestats(path1,path2):
    S1 = set(listfiles(path1))
    S2 = set(listfiles(path2))
    S = S1.intersection(S2)
    return S

def islarge(x1,x2):
    y1,q1 = x1
    y2,q2 = x2
    if y1 == 0:
        return False 
    elif y2 == 0 or q1-q2>2:
        return True
    else:
        return y1*10**(q1-q2)>y2

def summ(x1,x2):
    y1,q1 = x1
    y2,q2 = x2
    q = max(q1,q2)
    return (y1*10**(q1-q)+y2*10**(q2-q),q)

def runstest(al, threshold):
    N = len(al)
    n = 0
    for ch in al:
        if ch == 'm': 
            n += 1
    m = N - n   

    if min(m,n)<=threshold:
        return None,None
    R = 1
    for i in range(1,N):
        if (al[i-1] == "m" and al[i]!="m") or (al[i-1]!="m" and al[i]=="m"):
            R += 1
    #if R == 1:
    #    return None,None
    den = comb(N, n)
    q = int(math.log10(den))
    den = den/(10**q)
    for r in range(2,R+1):
        k = r//2
        if r%2==0:
            mk,nk = comb(m-1, k-1),comb(n-1, k-1)
            z = 2*mk*nk
        else:
            mk,nk,nkk,mkk = comb(m-1, k-1),comb(n-1, k-1),comb(n-1, k),comb(m-1, k)
            z = nkk*mk+mkk*nk 
        zq = int(math.log10(z))
        z/=10**zq
        y = z/den
        rq = zq-q
        if r == 2:
            prob = (y,rq)
        else:
            prob = summ(prob,(y,rq))
    y,q = prob
    logy = math.log10(y)
    x = round(logy - (0.5 if logy<0 else 0)) 
    y = y*10**(-x)
    q+=x
    return y,q

def countP(n, k):
     
    # Table to store results of subproblems
    dp = [[0 for i in range(k + 1)]
             for j in range(n + 1)]
 
    # Base cases
    for i in range(n + 1):
        dp[i][0] = 0
 
    for i in range(k + 1):
        dp[0][k] = 0
 
    # Fill rest of the entries in
    # dp[][] in bottom up manner
    for i in range(1, n + 1):
        for j in range(1, k + 1):
            if (j == 1 or i == j):
                dp[i][j] = 1
            else:
                dp[i][j] = (j * dp[i - 1][j] +
                                dp[i - 1][j - 1])
    
    return dp


def instest(al,threshold):
    Nm, Nd, Ns, ki = 0,0,0,0 #the number of matches, deletions, substitutions, insertions
    w = 0 #K is the number of non-empty bins
    atbin = True
    #need to count the number of matches + number of deletions + number of mismatches
    for x in al:
        if x == "-":
            ki+=1
            if atbin:
                atbin = False
                w+=1
        else:
            atbin = True
            if x == "m":
                Nm+=1
            elif x == "*":
                Ns+=1
            elif x == "/":
                Nd+=1
    N=Nm+Nd+Ns

    if ki<=threshold:
        return None
    if w == min(ki,N):
        return 1
    else:
        # Let's compute the probability of W<=w
        # We have to compute the stirling numbers
        dp = countP(ki,w) 
        ans = 0
        for W in range(1,w+1):
            logp=log(dp[ki][W]) + log(factorial(W)) + log(comb(N+1,W)) - ki*log(N+1)
            ans += exp(logp)
        return ans


def bfork(org,edt,k):
    n1 = len(org)
    n2 = len(edt)
    @lru_cache(None)
    def findb(i,j):
        if i>=n1 and j>=n2:
            return 0
        elif i>=n1:
            return ceil((n2-j)/k)
        elif j>=n2:
            return ceil((n1-i)/k)
        if org[i] == edt[j]:
            return findb(i+1,j+1)
        else:
            ans = float("inf")
            for k1 in range(k+1):
                for k2 in range(k+1):
                    if k1 == 0 and k2 == 0:
                        continue
                    ans = min(ans, findb(i+k1,j+k2))
        return ans+1
    return findb(0,0)


def smlestbfork(org,edt,k):
    dists = np.matrix(np.zeros((len(org)+1,len(edt)+1)))
    dists[0, 0] = 0
    for i in range(1, len(org)+1): dists[i, 0] = ceil(i/k)
    for i in range(1, len(edt)+1): dists[0, i] = ceil(i/k)
    for i in range(1, len(org)+1):
        for j in range(1, len(edt)+1):
            d = []
            org_i = org[0:i]
            edt_j = edt[0:j]
            for ii in range(0, min(k+1,i+1)):
                for jj in range(0, min(k+1,j+1)):
                    if ii != jj: d.append(dists[i-ii, j-jj] + 1)
                    elif ii != 0 or jj != 0:
                        if org_i[-ii:] == edt_j[-jj:]: d.append(dists[i-ii, j-jj])
                        else: d.append(dists[i-ii, j-jj] + 1)
            dists[i, j] = min(d)
    return dists[len(org),len(edt)]

def smlestbfork_dp(org,edt,k):
    dists = np.matrix(np.zeros((len(org)+1,len(edt)+1)))
    dists[0, 0] = 0
    for i in range(1, len(org)+1): dists[i, 0] = ceil(i/k)
    for i in range(1, len(edt)+1): dists[0, i] = ceil(i/k)
    for i in range(1, len(org)+1):
        for j in range(1, len(edt)+1):
            d = []
            if org[i-1] == edt[j-1]:
                dists[i,j] = dists[i-1,j-1]
            else:
                for ii in range(0, min(k+1,i+1)):
                    for jj in range(0, min(k+1,j+1)):
                        if ii==0 and jj==0:
                            continue 
                        d.append(dists[i-ii, j-jj] + 1)
                dists[i, j] = min(d)
    return int(dists[len(org),len(edt)])