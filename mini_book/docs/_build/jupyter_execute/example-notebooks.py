#!/usr/bin/env python
# coding: utf-8

# # Jupyter Notebook files
# 
# You can create content with Jupyter Notebooks. For example, the content for the current page is contained
# in {download}`this notebook file <./notebooks.ipynb>`. 
# 
# ```{margin}
# If you'd like to write in pure-text files, but still keep a notebook structure, you can write
# Jupyter Notebooks with MyST Markdown as well instead of `.ipynb`.
# See {doc}`myst-notebooks` for more details.
# ```
# 
# Jupyter Book supports all markdown that is supported by Jupyter Notebooks.
# This is mostly a flavor of markdown called
# [CommonMark Markdown](https://commonmark.org/) with minor modifications.
# For more information about writing Jupyter-flavored markdown in Jupyter Book,
# see {doc}`markdown`.
# 
# ## Code blocks and image outputs
# 
# Jupyter Book will also embed your code blocks and output in your book.
# For example, here's some sample Matplotlib code:

# In[1]:


from matplotlib import rcParams, cycler
import matplotlib.pyplot as plt
import numpy as np
plt.ion()


# In[2]:


# Fixing random state for reproducibility
np.random.seed(19680801)

N = 10
data = [np.logspace(0, 1, 100) + np.random.randn(100) + ii for ii in range(N)]
data = np.array(data).T
cmap = plt.cm.coolwarm
rcParams['axes.prop_cycle'] = cycler(color=cmap(np.linspace(0, 1, N)))


from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color=cmap(0.), lw=4),
                Line2D([0], [0], color=cmap(.5), lw=4),
                Line2D([0], [0], color=cmap(1.), lw=4)]

fig, ax = plt.subplots(figsize=(10, 5))
lines = ax.plot(data)
ax.legend(custom_lines, ['Cold', 'Medium', 'Hot']);


# Note that the image above is captured and displayed in your site.

# In[3]:


# Fixing random state for reproducibility
np.random.seed(19680801)

N = 10
data = [np.logspace(0, 1, 100) + .1*np.random.randn(100) + ii for ii in range(N)]
data = np.array(data).T
cmap = plt.cm.coolwarm
rcParams['axes.prop_cycle'] = cycler(color=cmap(np.linspace(0, 1, N)))


from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color=cmap(0.), lw=4),
                Line2D([0], [0], color=cmap(.5), lw=4),
                Line2D([0], [0], color=cmap(1.), lw=4)]

fig, ax = plt.subplots(figsize=(10, 5))
lines = ax.plot(data)
ax.legend(custom_lines, ['Cold', 'Medium', 'Hot'])
ax.set(title="Smoother linez");


# ```{margin} **You can also pop out content to the side!**
# For more information on how to do this,
# check out the {ref}`layout/sidebar` section.
# ```

# ## Removing content before publishing
# 
# You can also remove some content before publishing your book to the web. 
# For reference, {download}`you can download the notebook content for this page <notebooks.ipynb>`.

# In[4]:


thisvariable = "none of this should show up in the textbook"

fig, ax = plt.subplots()
x = np.random.randn(100)
y = np.random.randn(100)
ax.scatter(x, y, s=np.abs(x*100), c=x, cmap=plt.cm.coolwarm)
ax.text(0, .5, thisvariable, fontsize=20, transform=ax.transAxes)
ax.set_axis_off()


# You can **remove only the code** so that images and other output still show up.

# In[5]:


thisvariable = "this plot *will* show up in the textbook."

fig, ax = plt.subplots()
x = np.random.randn(100)
y = np.random.randn(100)
ax.scatter(x, y, s=np.abs(x*100), c=x, cmap=plt.cm.coolwarm)
ax.text(0, .5, thisvariable, fontsize=20, transform=ax.transAxes)
ax.set_axis_off()


# Which works well if you'd like to quickly display cell output without cluttering your content with code.
# This works for any cell output, like a Pandas DataFrame.

# In[6]:


import pandas as pd
pd.DataFrame([['hi', 'there'], ['this', 'is'], ['a', 'DataFrame']], columns=['Word A', 'Word B'])

