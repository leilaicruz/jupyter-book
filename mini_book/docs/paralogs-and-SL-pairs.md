---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: '0.8'
    jupytext_version: 1.4.2
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

```{code-cell} ipython3

import numpy as np
import scipy.io
import seaborn as sns
from scipy import stats, optimize, interpolate
import pandas as pd
from collections import defaultdict 
import math
import matplotlib.pyplot as plt
from scipy.stats import norm, lognorm
from scipy import stats
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import os, fnmatch
```

```{code-cell} ipython3
import os
script_dir = os.path.dirname('__file__') #<-- absolute dir the script is in
rel_path_SL = "datasets/data-synthetic-lethals.xlsx"

abs_file_path_domains = os.path.join(script_dir, rel_path_SL)

my_path_domains=abs_file_path_domains
data_domains=pd.read_excel(my_path_domains,header=0,index_col='Unnamed: 0')
data_sl=data_domains.dropna()

rel_path_query="datasets/paralogs_sL_from_query-genes.xlsx"
abs_file_path_query = os.path.join(script_dir, rel_path_query)

rel_path_target="datasets/paralogs_sL_from_target-genes.xlsx"
abs_file_path_target= os.path.join(script_dir, rel_path_target)

query_paralogs_pd=pd.read_excel(abs_file_path_query)
target_paralogs_pd=pd.read_excel(abs_file_path_target)
```

###  Build  a program that reads the paralogs of the que query gene and see if that paralog is also SL of the query gene , by inspecting if the paralog is present in the target genes of the SL database. 
- this is the first check to analyze if the reason why a SL pair shares domains is because they are also paralogs. 


```{code-cell} ipython3
# query_paralogs=pd.read_excel('datasets/paralogs_sL_from_query-genes.xlsx')
# target_paralogs=pd.read_excel('datasets/paralogs_sL_from_target-genes.xlsx')

query_paralogs_pd=query_paralogs.drop(columns='Unnamed: 0')
query_paralogs_pd.columns=['name-gene','name-paralogue']

target_paralogs_pd=target_paralogs.drop(columns='Unnamed: 0')
target_paralogs_pd.columns=['name-gene','name-paralogue']
```

```{code-cell} ipython3
query_paralogs_pd.head()
target_paralogs_pd.head()
```

```{code-cell} ipython3
indexes_sl_query=[]
for i in np.arange(0,len(query_paralogs_pd)):
    paralog_target=query_paralogs_pd[query_paralogs_pd['name-gene']==query_paralogs_pd['name-gene'][i]]['name-paralogue'].tolist()
    list_targets_sl=data_sl[data_sl['gene-query-name']==query_paralogs_pd['name-gene'][i]]['gene-target-name'].tolist()


    
    if paralog_target[0] in list_targets_sl:
        indexes_sl_query.append(query_paralogs_pd[query_paralogs_pd['name-paralogue']==paralog_target[0]].index[0])


indexes_sl_target=[]
for i in np.arange(0,len(target_paralogs_pd)): 
    
    paralog_query=target_paralogs_pd[target_paralogs_pd['name-gene']==target_paralogs_pd['name-gene'][i]]['name-paralogue'].tolist()
    list_queries_sl=data_sl[data_sl['gene-target-name']==target_paralogs_pd['name-gene'][i]]['gene-query-name'].tolist()


    if paralog_query[0] in list_queries_sl:
        indexes_sl_target.append(target_paralogs_pd[target_paralogs_pd['name-paralogue']==paralog_query[0]].index[0])

```

### Putting 1's if the paralog pair is also SL

```{code-cell} ipython3
sL_values=np.zeros_like(query_paralogs_pd['name-gene'])
for i in np.arange(0,len(query_paralogs_pd)):
    if i in indexes_sl_query:
        sL_values[i]=1
query_paralogs_pd['sL']=sL_values

sL_values=np.zeros_like(target_paralogs_pd['name-gene'])
for i in np.arange(0,len(target_paralogs_pd)):
    if i in indexes_sl_target:
        sL_values[i]=1
target_paralogs_pd['sL']=sL_values
```

```{code-cell} ipython3
paralogs_sl_pd=pd.concat([query_paralogs_pd,target_paralogs_pd],axis=0)
```

```{code-cell} ipython3
sl_that_are_paralogs=paralogs_sl_pd[paralogs_sl_pd['sL']==1]
```

```{code-cell} ipython3
sl_that_are_paralogs.set_index(np.arange(0,len(sl_that_are_paralogs)))
```

```{code-cell} ipython3
sl_that_are_paralogs.iloc[0,0:2].tolist()
```

```{code-cell} ipython3
a= ['BUD6', 'CDC42']
b= ['CDC42', 'BUD6']
set(a)==set(b)
```

```{code-cell} ipython3
## use this to compare with the pairs of SL that have shared domains from FEATUREPROCESSINg script  
```

```{code-cell} ipython3
pairs_sL=np.load('../pairs-sL-that-share-domains.npy')
```

```{code-cell} ipython3
len(pairs_sL)
```

```{code-cell} ipython3
len(sl_that_are_paralogs)
```

```{code-cell} ipython3
shared_sL_paralogs=[]
for i in np.arange(0,len(sl_that_are_paralogs)):
    for j in np.arange(0,len(pairs_sL)):
        if set(sl_that_are_paralogs.iloc[i,0:2].tolist())==set(pairs_sL[j]):
            shared_sL_paralogs.append(pairs_sL[j])

```

```{code-cell} ipython3
:tags: []

print('The contribution of paralogs to the SL pairs that shared domains is =', 100*len(shared_sL_paralogs)/len(pairs_sL),'%')

print('The contribution of paralogs to the total number of SL pairs is  =', 100*len(sl_that_are_paralogs)/17871,'%')

print('The number of SL that share domains out of the total number of SL pairs is =',100*len(pairs_sL)/17871,'%')
```

```{code-cell} ipython3

```
