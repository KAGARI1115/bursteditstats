#%%
from func import islarge,runstest,instest
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt 
#%%

#apply the runs test, insertion test and deletion test on the genome dataset,
#and plot the histgram for p-values

diry = "genome/guppy/"

p_runs_al,p_runs_file, p_runs_al_list = (0,0),"",[]
p_ins,p_ins_file, p_ins_list = 0,"",[]
p_deltest,p_del_file, p_deltest_list = (0,0),"",[]

allim, inslim, dellim = 20,20,20

alcnt, alpass = 0,0
inscnt,inspass = 0,0
delcnt,delpass=0,0

threshold = 0
ifplot = True
for f in range(1000):
    filename = diry+"j"+str(f)+".txt" 
    al_file = open(filename , "r")
    als = al_file.readlines()
    al_file.close()

    org = als[-1]
    al = als[-2][:-1]
    edt = als[-3][:-1]

    true_al = ""
    del_al = ""
    for i,x in enumerate(al):
        if x == "|":
            true_al = true_al+"m"
            del_al=del_al+"m"
        elif x==".":
            true_al = true_al+"*"
            del_al=del_al+"*"
        elif x=="-":
            if org[i]=="-":
                true_al = true_al+"-"
            else:
                true_al=true_al+"/"
                del_al=del_al+"/"

    # apply runs-test on true al
    (y,p) = runstest(true_al, 2*threshold)
    if y:
        p_runs_al_list.append(y*10**(p))
        if y*10**(p)<=allim/100:
            alcnt+=1
        if y*10**(p)<=0.05:
            alpass+=1
        if islarge((y,p),p_runs_al):
            p_runs_al = (y,p)
            p_runs_file=filename
    
    # apply ins-test on true al
    pins = instest(true_al,5)
    if pins:
        p_ins_list.append(pins)
        if pins<=inslim/100:
            inscnt+=1
        if pins<=0.05:
            inspass+=1
        if pins>p_ins:
            p_ins = pins
            p_ins_file = filename

    #apply runs-test on del and subs
    (y,p) = runstest(del_al, threshold)
    if y:
        p_deltest_list.append(y*10**(p))
        if y*10**(p)<=dellim/100:
            delcnt+=1
        if y*10**(p)<=0.05:
            delpass+=1
        if islarge((y,p),p_deltest):
            p_deltest = (y,p)
            p_del_file=filename

print(p_runs_al,p_runs_file,len(p_runs_al_list) ,alcnt,alpass, alpass/len(p_runs_al_list))
print(p_ins,p_ins_file,len(p_ins_list) ,inscnt,inspass,inspass/len(p_ins_list))
print(p_deltest,p_del_file,len(p_deltest_list) ,delcnt,delpass,delpass/len(p_deltest_list))
#%%
if ifplot:
    plt.figure(figsize=(7,6))
    hist_al = sns.histplot(data=np.array(p_runs_al_list),bins=[x/allim for x in range(allim+1)] )
    plt.xlabel("$p$-values of runs-test on alignment",fontsize = 20)
    plt.ylabel('count', fontsize = 20)
    plt.xticks([x/allim for x in range(allim+1) if x%2==0],fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.savefig("genome/p_runs_al_genome.png")

    plt.figure(figsize=(7,6))
    hist_ins = sns.histplot(data=np.array(p_ins_list),bins=[x/inslim for x in range(inslim+1)] )
    plt.xlabel("$p$-values of ins-test",fontsize = 20)
    plt.ylabel('count', fontsize = 20)
    plt.xticks([x/inslim for x in range(inslim+1) if x%2==0],fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.savefig("genome/p_ins_genome.png")

    plt.figure(figsize=(7,6))
    hist_del = sns.histplot(data=np.array(p_deltest_list),bins=[x/dellim for x in range(dellim+1)] )
    plt.xlabel("$p$-values of subdel-test",fontsize = 20)
    plt.ylabel('count', fontsize = 20)
    plt.xticks([x/dellim for x in range(dellim+1) if x%2==0],fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.savefig("genome/p_deltest_genome.png")



    # y,q = runstest(al,0)
    # print("filename:",file_name)
    # print("one-sided p-value=",y,"e",q)
    # print("\n")


#%%
# filename = ""
# mp_value = (0,0)
# mfilename = ""
# for i,al in enumerate(als):
#     if i % 2 == 0:
#         filename = al
#         print(filename)
#     else:
#         y,q = runstest(al,10)
#         if y:
#             if islarge((y,q),mp_value):
#                 mp_value = (y,q)
#                 mfilename = filename
# print(mp_value)
# print(mfilename)