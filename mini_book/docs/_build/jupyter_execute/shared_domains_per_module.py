#!/usr/bin/env python
# coding: utf-8

# # What are the most shared domains among proteins that belong to a certain biological process module? 
# 
#  In this notebook , I inspect the genes that belong to certain GO terms , search for their domains , and sort by the most frequent domain to the lowest.

# In[1]:


from intermine.webservice import Service
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


# In[2]:


import os
script_dir = os.path.dirname('__file__') #<-- absolute dir the script is in
rel_path_domains="datasets/proteins-domains-from-Pfam.xlsx"


abs_file_path_domains = os.path.join(script_dir, rel_path_domains)

# os.chdir('../') #<-- for binder os.chdir('../')
# os.chdir('mini_book/docs/')

my_path_domains=abs_file_path_domains
data_domains=pd.read_excel(my_path_domains,header=0,index_col='Unnamed: 0')
data_domains=data_domains.dropna()

