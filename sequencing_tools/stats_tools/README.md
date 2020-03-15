# sequencing_tools.stats_tools #


## multiple testing p-value ##

Benjamini-Hochberg p-value correction for multiple hypothesis testing.
adapted from https://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python/21739593

Parameter:

- pvalue_array: list of pvalue

Return:

- list of adjusted pvalue 

Usage:

```
from sequencing_tools.stats_tools import p_adjust
padj = p_adjust(pvalue_array)
```


## Vectorized binomial test ##

Calculate p-values for multiple binomial tests, let's say you have 5 coins and each of them were flipped 10 times, and you want to test if they are unbiased coin

input:
    failed_test: integer array indicating number times you get "tails" for each coin
    total_test: integer array indicating number of total flips you have conducted for each coin
    expected_p: the expected p-value to be test against

return:
    list of p-values
ectorized binomial test:


usage:
```
from sequencing_tools.stats_tools import binom_test
failed_test = [0,5,6,19,1]
total_test = [10,6,7,20,10]
pvalue_array = binom_test(failed_test, total_test, expected_p = 0.5)
```


## Mean ##

Fast numerical mean calculator, works with number generator too

usage: cy_mean(list_of_numbers)
    return: mean of the list

Usage:
```
from sequencing_tools.stats_tools import cy_mean
a = range(10)
cy_mean(a)
```

Test:

    In [1]: from sequencing_tools.stats_tools import cy_mean
    In [2]: import numpy as np
    In [3]: a = range(10)
    In [4]: b = np.array(a)
    In [5]: %timeit b.mean()
    6.37 µs ± 323 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)
    In [6]: %timeit cy_mean(a)
    385 ns ± 7.75 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)


## levenshtein distance ##

Calculating Levenshtein distance from two strings
algorithm from: http://rosettacode.org/wiki/Levenshtein_distance#Python

input:
    string1
    string2

return:
    edit distance: the edit distance between two string

Usage:
```
from sequencing_tools.stats_tools import levenshtein_distance
string1 = 'AACGGACT'
string2 = 'ACTGACGG'
levenshtein_distance(string1, string2) #4
```


## hamming distance ##

Calculating hamming distance from two strings

input:
    string1: must be same length as string2
    string2: must be same length as string1

return:
    edit distance: the edit distance between two string

Usage:
```
from sequencing_tools.stats_tools import hamming_distance
string1 = 'AACGGACT'
string2 = 'ACTGACGG'
hamming_distance(string1, string2) #6
```

## Normalize count ##
Normalize a RNA-seq count matrix using DESeq2 method. The DESeq2 size factor algorithm is described [here](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106#Sec22).


input:
    pandas dataframe with shape(m,n): m is gene count, n is sample number
    return_sf: return size factors (n)

output:
    np.ndarray(m,n)

```
df
```
| Name           |   sample1 |   sample2 |   sample3 |   sample4 |
|:---------------|------------:|------------:|------------:|------------:|
| hsa-miR-16-5p  |      242711 |      285618 |      457227 |      363398 |
| hsa-miR-486-5p |      193648 |      238126 |      316513 |      376722 |
| hsa-let-7b-5p  |       83115 |       84112 |      116224 |      126184 |
| hsa-miR-223-3p |       81857 |       56294 |       31904 |       61472 |
| hsa-miR-92a-3p |       54848 |       50318 |       56969 |       65637 |
| hsa-miR-25-3p  |       20500 |       39468 |       35198 |       49625 |
| hsa-let-7a-5p  |       46059 |       39187 |       68195 |       54559 |
| hsa-miR-451a   |       16495 |       35044 |       18004 |       23891 |
| hsa-miR-191-5p |       17823 |       23148 |       13376 |       21004 |
| hsa-miR-30d-5p |       17774 |       22308 |       14948 |       24809 |

```
from sequencing_tools.stats_tools import normalize_count
norm_df, sf = normalize_count(df.set_index('Name'), return_sf = True)
norm_df
```
| Name           |  sample1 |   sample2 |   sample3 |   sample4  |
|:---------------|------------:|------------:|------------:|------------:|
| hsa-miR-16-5p  |    280758   |    299162   |    451377   |    316506   |
| hsa-miR-486-5p |    224004   |    249418   |    312464   |    328110   |
| hsa-let-7b-5p  |     96143.9 |     88100.7 |    114737   |    109901   |
| hsa-miR-223-3p |     94688.7 |     58963.5 |     31495.8 |     53539.7 |
| hsa-miR-92a-3p |     63445.8 |     52704.1 |     56240.1 |     57167.3 |
| hsa-miR-25-3p  |     23713.5 |     41339.6 |     34747.7 |     43221.5 |
| hsa-let-7a-5p  |     53279.1 |     41045.3 |     67322.5 |     47518.8 |
| hsa-miR-451a   |     19080.7 |     36705.8 |     17773.7 |     20808.1 |
| hsa-miR-191-5p |     20616.9 |     24245.7 |     13204.9 |     18293.7 |
| hsa-miR-30d-5p |     20560.2 |     23365.9 |     14756.8 |     21607.7 |

```
sf
array([0.86448575, 0.95472601, 1.01295971, 1.14815656])
```


