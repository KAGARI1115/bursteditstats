#%%
from collections import defaultdict
from func import bfork,filestats,smlestbfork,align,smlestbfork_dp
import pandas as pd
import sys
sys.setrecursionlimit(100000)

#%%
# f = open('bk.csv','r')
# filenames = set([x.split(',')[0] for x in f.readlines()])


#%%
mypath1 = "./bash-5.0"
mypath2 = "./bash-5.1"

c=0
length = 0
K = 10

for f in filestats(mypath1, mypath2):
    filename1 = mypath1 +f
    filename2 = mypath2 +f
    a = list(open(filename1,"rb").read())
    b = list(open(filename2,"rb").read())
    if a!=b and len(a)<3000:#1000<=len(a)<1500 and f not in filenames:
        print("filename", f)
        print(len(a))
        dict = {'name':f,'file size':len(a)}
        for k in range(1,K+1):
            print(k)
            dict['k='+str(k)]=smlestbfork_dp(a,b,k)
        row = pd.DataFrame(dict, index = [f])
        if c==0:
            row.to_csv('bk_new.csv',mode='a',index=False)
        else:
            row.to_csv('bk_new.csv',mode='a',index=False, header=False)
        c+=1
        
