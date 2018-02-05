#!/usr/bin/env python

import numpy as np
cimport numpy as np
from cython cimport floating


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