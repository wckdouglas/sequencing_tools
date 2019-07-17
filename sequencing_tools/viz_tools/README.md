
# sequencing_tools.viz_tools #


## Upset graph ##


```python
%matplotlib inline
%load_ext autoreload
%autoreload 2

import matplotlib.pyplot as plt
from sequencing_tools.viz_tools import plot_upset
import numpy as np
from itertools import combinations
import pandas as pd
```

    The autoreload extension is already loaded. To reload it, use:
      %reload_ext autoreload


Let's simulate a dataset:


```python
var = 'ABCD'
encoder = {e:i for i, e in enumerate(var)}
matrix = []
for comb1, comb2 in combinations(var, 2):
    row = np.zeros(len(var))
    row[encoder[comb1]] += 1
    row[encoder[comb2]] += 1
    matrix.append(row)
upset_df = pd.DataFrame(matrix,columns=list(var)) \
    .assign(count = lambda d: np.random.random_integers(0,high = 10,size=d.shape[0]))\
    .reset_index()
upset_df
```

    /home/wckdouglas/miniconda3/lib/python3.6/site-packages/ipykernel_launcher.py:10: DeprecationWarning: This function is deprecated. Please call randint(0, 10 + 1) instead
      # Remove the CWD from sys.path while we load stuff.





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>index</th>
      <th>A</th>
      <th>B</th>
      <th>C</th>
      <th>D</th>
      <th>count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>5</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>10</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>7</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>8</td>
    </tr>
    <tr>
      <th>4</th>
      <td>4</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>2</td>
    </tr>
    <tr>
      <th>5</th>
      <td>5</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>6</td>
    </tr>
  </tbody>
</table>
</div>



The ```plot_upset``` function requires a input of ```plt.figure``` as placeholder and a pandas dataframe as upset_df


```python
fig = plt.figure()
plot_upset(fig, upset_df.sort_values('count'), ylab='count', matrix_to_plot_ratio=0.4, fontsize=20)
```


![png]((https://raw.githubusercontent.com/wckdouglas/sequencing_tools/master/img/upset.png)


## Correlation matrix scatter plot ##


```python
from sequencing_tools.viz_tools import cor_plot
import seaborn as sns
sns.set_style('white')

n_sample = 5
d = np.random.rand(10,n_sample)
d = pd.DataFrame(d, columns=['Sample %i' %i for i in range(n_sample)])
fig = plt.figure(figsize=(8,8))
cor_plot(d, fig)
```


![png]((https://raw.githubusercontent.com/wckdouglas/sequencing_tools/master/img/scatter.png)

