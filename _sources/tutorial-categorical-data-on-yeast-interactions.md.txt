---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

(tutorial-categorical-data-on-yeast-interactions)=

# Tutorial on how to handle categorical data from [here](https://www.datacamp.com/community/tutorials/categorical-data)


```{code-cell} python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict 
import seaborn as sns
```


```{code-cell} python3
data=pd.read_excel(r'C:\Users\linigodelacruz\Documents\PhD_2018\Documentation\Calculations\data_BioGrid\data-BioGrid-Yeast.xlsx',header=0)

```


```{code-cell} python3
data = data[['gene-query-name','gene-target-name','interaction-type', 'paper-source']]
```


```{code-cell} python3
data.head()
```


```{code-cell} python3
print(data.info())
```

   

As you will only be dealing with categorical features in this tutorial, it's better to filter them out. You can create a separate DataFrame consisting of only these features by running the following command. The method ```.copy()``` is used here so that any changes made in new DataFrame don't get reflected in the original one.


```{code-cell} python3
cat_data = data.select_dtypes(include=['object']).copy()
```


```{code-cell} python3
cat_data.head()
```




One of the most common data pre-processing steps is to check for null values in the dataset. You can get **the total number of missing values in the DataFrame** by the following one liner code:


```{code-cell} python3
if cat_data.isnull().values.sum()==0:
    print('Hooray!! There are no NaN values in the dataframe')
else:
    print(cat_data.isnull().values.sum())
```

    Hooray!! There are no NaN values in the dataframe
    

Let's also check the column-wise distribution of null values:


```{code-cell} python3
print(cat_data.isnull().sum())
```

Another **Exploratory Data Analysis (EDA)** step that you might want to do on categorical features is the frequency distribution of categories within the feature, which can be done with the ```.value_counts()``` method as described earlier.


```{code-cell} python3
print(cat_data['interaction-type'].value_counts())
```

     


```{code-cell} python3
print(cat_data['gene-query-name'].value_counts())
```

      

To know **the count of distinct categories within the feature** you can chain the previous code with the ```.count()``` method:


```{code-cell} python3
print(cat_data['interaction-type'].value_counts().count(),cat_data['gene-query-name'].value_counts().count(),cat_data['gene-target-name'].value_counts().count())
```

    

Below is a basic template to plot a barplot of the frequency distribution of a categorical feature using the seaborn package, which shows the frequency distribution of the carrier column. You can play with different arguments to change the look of the plot.


```{code-cell} python3
carrier_count = cat_data['interaction-type'].value_counts()
sns.set(style="darkgrid")
sns.barplot(carrier_count.index, carrier_count.values, alpha=0.9)
plt.title('Frequency Distribution of interaction types in BioGrid')
plt.ylabel('Number of Occurrences', fontsize=12)
plt.xlabel('Interaction types', fontsize=12)
plt.xticks(rotation=30)
```




## Encoding Categorical Data

The techniques that you'll cover are the following:

1. Replacing values
2. Encoding labels
3. One-Hot encoding
4. Binary encoding

### Replace Values

Let's start with the most basic method, which is just replacing the categories with the desired numbers. This can be achieved with the help of the ```replace()``` function in pandas. The idea is that you have the liberty to choose whatever numbers you want to assign to the categories according to the business use case.

You will store the category names in a list called ```labels``` and then ```zip``` it to a sequence of numbers and iterate over it. The final dictionary will organize the labels in alphabetical order. 


```{code-cell} python3
labels = cat_data['interaction-type'].astype('category').cat.categories.tolist()
replace_map_comp = {'interaction-type' : {k: v for k,v in zip(labels,list(range(1,len(labels)+1)))}}
```


```{code-cell} python3
print(replace_map_comp)
```

    

Throughout this tutorial, you will be making a copy of the dataset via the ```.copy()``` method to practice each encoding technique to ensure that the original DataFrame stays intact and whatever changes you are doing happen only in the copied one.


```{code-cell} python3
cat_data_replace = cat_data.copy()
```


```{code-cell} python3
cat_data_replace.replace(replace_map_comp, inplace=True)

print(cat_data_replace.head())
```
    

As you can observe, you have encoded the categories with the mapped numbers in your DataFrame.

## Label encoding

Another approach is to encode categorical values with a technique called "label encoding", which allows you to convert each value in a column to a number. Numerical labels are always between 0 and n_categories-1

You can do label encoding via attributes ```.cat.codes``` on your DataFrame's column.


```{code-cell} python3
cat_data_lc=cat_data.copy().astype('category') # In general converting to a  category variable is much faster and handy that leaves them as object

```


```{code-cell} python3
cat_data_lc['interaction-type'] = cat_data_lc['interaction-type'].cat.codes
```


```{code-cell} python3
cat_data_lc.head() #alphabetically labeled from 0 to number of categories : 27
```




Sometimes, **you might just want to encode a bunch of categories within a feature to some numeric value and encode all the other categories to some other numeric value**.

You could do this by using numpy's ```where()``` function like shown below. 
Example: You will encode all the synthetic lethals to value 1 and other types to value 0. This will create a new column in your DataFrame with the encodings. Later, if you want to drop the original column, you can do so by using the ```drop()``` function in pandas.


```{code-cell} python3
cat_data_specific = cat_data.copy()
cat_data_specific['SL-code'] = np.where(cat_data_specific['interaction-type'].str.contains('Lethality'), 1, 0)

cat_data_specific.head()
```



```{code-cell} python3
cat_data_specific[cat_data_specific['SL-code']==1].head()
```




## You can achieve the same label encoding using `scikit-learn's, LabelEncoder` 



```{code-cell} python3
cat_data_sklearn = cat_data.copy()

from sklearn.preprocessing import LabelEncoder

lb_make = LabelEncoder()
cat_data_sklearn['type_code'] = lb_make.fit_transform(cat_data['interaction-type'])

cat_data_sklearn.tail() #Results in appending a new column to df
```




Label encoding is pretty much intuitive and straight-forward and may give you a good performance from your learning algorithm, but it has as disadvantage that the numerical values can be misinterpreted by the algorithm. 

To solve this issue there is another popular way to encode the categories via something called one-hot encoding.


## One-Hot encoding

The basic strategy is to convert each category value into a new column and assign a 1 or 0 (True/False) value to the column. This has the benefit of not weighting a value improperly.

There are many libraries out there that support one-hot encoding but the simplest one is using pandas' ```.get_dummies()``` method.

This function is named this way because it creates dummy/indicator variables (1 or 0). There are mainly three arguments important here, the first one is the DataFrame you want to encode on, second being the columns argument which lets you specify the columns you want to do encoding on, and third, the prefix argument which lets you specify the prefix for the new columns that will be created after encoding.


```{code-cell} python3
cat_data_onehot = cat_data.copy()
cat_data_onehot = pd.get_dummies(cat_data_onehot, columns=['interaction-type'])


cat_data_onehot.head()
```



```scikit-learn``` also supports one hot encoding via ```LabelBinarizer``` and ```OneHotEncoder``` in its preprocessing module (check out the details here). Just for the sake of practicing you will do the same encoding via ```LabelBinarizer```:


```{code-cell} python3
cat_data_onehot_sklearn = cat_data.copy()

from sklearn.preprocessing import LabelBinarizer

lb = LabelBinarizer()
lb_results = lb.fit_transform(cat_data_onehot_sklearn['interaction-type'])
lb_results_df = pd.DataFrame(lb_results, columns=lb.classes_)

lb_results_df.head()
```







```{code-cell} python3
result_df = pd.concat([cat_data_onehot_sklearn, lb_results_df], axis=1)

result_df.head()
```


**While one-hot encoding solves the problem of unequal weights given to categories within a feature**, it is not very useful when there are many categories, as that will result in formation of as many new columns, which can result in the curse of dimensionality. The concept of the “curse of dimensionality” discusses that in high-dimensional spaces some things just stop working properly.






