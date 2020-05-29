(machine-learning-replication-protein-domains)=

# Example

- Replication of results from paper: "Predicting yeast synthetic lethal genetic interactions using protein domains"
    - Authors: Bo Li, Feng Luo,School of Computing,Clemson University,Clemson, SC, USA
    - e-mail: bol, luofeng@clemson.edu
    - year:2009

%matplotlib inline
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict 
import seaborn as sns
import matplotlib.cm as cm
import scipy as scipy
import random