#%%
import pandas as pd 
from math import log2 
import numpy as np

df = pd.read_excel('answer_2.xlsx')

print(df)
#%%
R = pd.DataFrame()
R['file name'] = df['file name']



q = 4
for k in range(1,10):
    col = "k = "+str(k)
    R[col] = df[col].apply(lambda x: 2*x*log2(200+x*k)+4*x*log2(k+1)+2*k*x*log2(q)).round(2)

# R['min'] = R.min(axis = 1)
# usableR = R[R['min']<400]
# usabledf = df[R['min']<400]

# # Medians = usabledf.median()
# # print(Medians)
# # half = [2*x*log2(200+x*(k+1))+4*x*log2(k+2)+2*(k+1)*x*log2(q) for k,x in enumerate(Medians.values)]
# # print(half)

# for i in range(1,21):
#     a = i*0.05
#     a = round(a,2)
#     #quantiles = df.quantile(q=a, interpolation='higher')
#     quantiles = usabledf.quantile(q=a, interpolation='higher')
#     qut = [2*x*log2(200+x*(k+1))+4*x*log2(k+2)+2*(k+1)*x*log2(q) for k,x in enumerate(quantiles.values[0:1])]
#     mink = np.argmin(qut)+1
#     print('percent = ',a) 
#     print('k =',mink, 'b = ', quantiles.values[mink-1], 'redundancy = ',round(qut[mink-1],2))
#     print('\n')
