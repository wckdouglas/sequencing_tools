from __future__ import print_function

import logging

import numpy as np
import pandas as pd
from scipy.stats import binom

from ..utils import SeqUtilsError

cimport numpy as np
logging.basicConfig(level = logging.INFO)
logger = logging.getLogger('Stat tool')

cpdef np.ndarray p_adjust(pvalue):
    '''
    Benjamini-Hochberg p-value correction for multiple hypothesis testing.
    adapted from https://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python/21739593

    Usage::

        padj = p_adjust(pvalue_array)

    Args:
        pvalue_array: list of pvalue

    Returns:
        list: list of adjusted pvalue 

    An equivalent in R::

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
    if pvalue.ndim != 1:
        raise SeqUtilsError('Only accept 1D array')

    descending_p_order = pvalue.argsort()[::-1]
    in_order = descending_p_order.argsort()
    steps = float(len(pvalue)) / np.arange(len(pvalue), 0, -1)
    adjusted_p = np.minimum(1.0, np.minimum.accumulate(steps * pvalue[descending_p_order]))
    logger.info('Corrected %i p-values' %len(pvalue))
    return adjusted_p[in_order]

def binom_test(success_test, total_test, expected_p = 0.5):
    '''
    Vectorized binomial test

    Usage::

        success_test = [10,10,10]
        total_test = [100,100,100]
        binom_test(success_test, total_test, expected_p = 0.5)
        #array([1.53164509e-17, 1.53164509e-17, 1.53164509e-17])

    Args:
        success_test: number of success
        total_test: number of total trials
        expected_p: expected probability of success

    Returns:
        list: list of p-values
        
    '''
    if len(success_test)!=len(total_test):
        raise SeqUtilsError('Wrong length of vector!')
    ps = binom.cdf(success_test,n=total_test,p=expected_p)
    return ps

cpdef double cy_mean(xs):
    '''
    Fast numerical mean calculator, works with number generator too

    Example::

        from sequencing_tools.stats_tools import cy_mean
        import numpy as np
        a = range(10)
        b = np.array(a)
        %timeit b.mean()
        # 6.37 µs ± 323 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)
        %timeit cy_mean(a)
        # 385 ns ± 7.75 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)

    Usage:: 
        
        cy_mean(list_of_numbers)

    Args:
        xs: list of numbers

    Returns: 
        float: Mean of the list
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

    usage:: 
        
        levenshtein_distance(string1, string2)
    
    Args:
        string1: first string in the comparison
        string2: second string in the comparison

    Returns:
        int: the edit distance between two string
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

    Calculating hamming distance from two strings, the two strings has to be in the same length

    Usage:: 
        
        hamming_distance(string1, string2)
    
    Args:
        string1: first string in the comparison
        string2: second string in the comparison

    Returns:
        int: the hamming edit distance between two strin

    '''

    cdef:
        str i, j
        int hamming = 0
    if len(s1) != len(s2):
        raise SeqUtilsError( 'Wrong barcode extraction')

    for i, j in zip(s1, s2):
        if i != j:
            hamming += 1

    return hamming


def normalize_count(count_mat, return_sf = False):
    '''
    `DESeq2 size factor <https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106#Sec22>`_

    Args:
        count_mat: pandas dataframe with shape(m,n): m is gene count, n is sample number
        return_sf: boolean, if it needs returning size factors 

    Returns:
        np.ndarray(m,n): normalized count matrix
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
    logger.info('Calculated size factor for %i samples' %len(sf))
    norm_count = count_mat/sf
    if return_sf:
        return norm_count, sf
    else:
        return norm_count


class Bootstrap:
    """
    Args:
        seed: seed for random number generater
    """
    def __init__(self, seed=123):
        self.rng = np.random.RandomState(seed) 

    def generate(self, xs, group_size=100, n_boots = 100):
        '''
        boostrap 1d array

        Usage::

            xs = np.arange(100)
            bs = Bootstrap(seed=123)
            for idx in bs.generate(xs, group_size=50, n_boots=10):
                print(xs[idx].mean())

        Args:
            xs: 1d np.array
            group_size: number of values in each bootstrap iteration
            n_boots: how many bootstrap groups
        Returns:
            iterator: bootstrapped
        '''
        xs = np.array(xs)
        total_size = xs.shape[0]
        logger.info('Total size for bootstrap: %i' %total_size)
        if group_size > total_size:
            #raise SeqUtilsError('Group size > input array size')
            raise SeqUtilsError('Group size > input array size')
    
        for i in range(n_boots):
            idx = self.rng.randint(0, total_size, group_size)
            yield idx
