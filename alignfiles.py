#%%
from collections import defaultdict
from func import align,filestats
#%%

dedits = defaultdict(int)

def editsfar(align):
    last = 0
    for i,x in enumerate(align):
        if x != "m":
            dedits[i-last-1]+=1
            last = i 

#%%
mypath1 = "./bash-5.0"
mypath2 = "./bash-5.1"
c=0
length = 0
#fw = open('./alignment.txt','w')
for f in filestats(mypath1, mypath2):
    filename1 = mypath1 +f
    filename2 = mypath2 +f
    a = list(open(filename1,"rb").read())
    b = list(open(filename2,"rb").read())
    if a!=b and len(a)<20000:
        c+=1
        # print("filename", f)
        # print(len(a))
        # al, _ = align(a,b)
        # fw.write(f)
        # fw.write("\n")
        # fw.write(al)
        # fw.write("\n")
print(c)
