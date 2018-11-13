from __future__ import print_function
import numpy as np
cimport numpy as np
from cython cimport floating
from scipy.stats import binom
import pandas as pd

cpdef np.ndarray p_adjust(pvalue):
    '''
    Benjamini-Hochberg p-value correction for multiple hypothesis testing.
    adapted from https://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python/21739593

    Usage:

    padj = p_adjust(pvalue_array)

    Parameter:

    - pvalue_array: list of pvalue

    Return:

    - list of adjusted pvalue 
    '''

    '''    
    in R:

        nm <- names(p)
        vp <- as.numeric(p)
        p0 <- setNames(p, nm)
        lp <- length(p)
        
        i <- lp:1L
        o <- order(p, decreasing = TRUE)
        ro <- order(o)
        pmin(1, cummin(n/i * p[o]))[ro]
    '''

    cdef:
        np.ndarray  steps, adjusted_p
        np.ndarray[long, ndim=1] descending_p_order, in_order

    pvalue = np.array(pvalue, dtype='float')
    pvalue[np.isnan(pvalue)] = 1
    assert pvalue.ndim == 1, 'Only accept 1D array'
    descending_p_order = pvalue.argsort()[::-1]
    in_order = descending_p_order.argsort()
    steps = float(len(pvalue)) / np.arange(len(pvalue), 0, -1)
    adjusted_p = np.minimum(1.0, np.minimum.accumulate(steps * pvalue[descending_p_order]))
    return adjusted_p[in_order]

def binom_test(success_test, total_test, expected_p = 0.5):
    '''
    Vectorized binomial test:

    usage:

        binom_test(failed_test, total_test, expected_p = 0.5)
    
    return:
        list of p-values
        
    '''
    assert len(success_test)==len(total_test), 'Wrong length of vector!'
    ps = binom.cdf(success_test,n=total_test,p=expected_p)
    return ps

cpdef double cy_mean(xs):
    '''
    Fast numerical mean calculator, works with number generator too


    In [1]: from sequencing_tools.stats_tools import cy_mean
    In [2]: import numpy as np
    In [3]: a = range(10)
    In [4]: b = np.array(a)
    In [5]: %timeit b.mean()
    6.37 µs ± 323 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)
    In [6]: %timeit cy_mean(a)
    385 ns ± 7.75 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)

    usage: cy_mean(list_of_numbers)
    return: mean of the list
    '''
    cdef:
        double x
        int counter = 0
        double sum_x
        double res
    
    for x in xs:
        counter += 1
        sum_x += x

    res = sum_x / counter
    return res

       
cpdef int levenshtein_distance(str s1, str s2):
    '''
    Calculating Levenshtein distance from two strings
    algorithm from: http://rosettacode.org/wiki/Levenshtein_distance#Python

    usage: levenshtein_distance(string1, string2)
    ==============================
    Parameter:

    string1
    string2

    return:
    edit distance: the edit distance between two string

    '''

    cdef:
        list distance, new_distances
        int index1, index2
        str char1, char2
        int len_s1 = len(s1)
        int len_s2 = len(s2)


    if len_s1 > len_s2:
        s1, s2 = s2, s1
    distances = range(len_s1 + 1)
    for index2, char2 in enumerate(s2):
        new_distances = [index2+1]
        for index1, char1 in enumerate(s1):
            if char1 == char2:
                new_distances.append(distances[index1])
            else:
                new_distances.append(1 + min((distances[index1],
                                             distances[index1+1],
                                             new_distances[-1])))
        distances = new_distances
    return distances[-1]


cpdef int hamming_distance(str s1, str s2):
    '''
    Calculating hamming distance from two strings

    usage: hamming_distance(string1, string2)
    ==============================
    Parameter:

    string1
    string2

    has to be same length

    return:
    edit distance: the edit distance between two string
    ===============================
    '''

    cdef:
        str i, j
        int hamming = 0
    assert len(s1) == len(s2), 'Wrong barcode extraction'

    for i, j in zip(s1, s2):
        if i != j:
            hamming += 1

    return hamming


def normalize_count(count_mat, return_sf = False):
    '''
    DESeq2 size factor:
        https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106#Sec22

    input:
        pandas dataframe with shape(m,n): m is gene count, n is sample number
        return_sf: return size factors (n)

    output:
        np.ndarray(m,n)
    '''

    cols = count_mat.columns.tolist()
    log_row_geomean = count_mat.transform(np.log)\
            .mean(axis=1, skipna=True)
    finite_mean = np.isfinite(log_row_geomean)
    
    
    def get_size_factor(sample_cnts):
        ratio = np.log(sample_cnts) - log_row_geomean
        ratio = ratio[(finite_mean) & (sample_cnts>0)]
        ratio = np.median(ratio)
        return ratio
    
    sf = np.apply_along_axis(get_size_factor, 0, count_mat)
    sf = np.exp(sf)
    norm_count = count_mat/sf
    if return_sf:
        return norm_count, sf
    else:
        return norm_count
