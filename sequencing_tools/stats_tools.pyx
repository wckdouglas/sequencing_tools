#!/usr/bin/env python

import numpy as np
cimport numpy as np


cpdef np.ndarray p_adjust(np.ndarray pvalue):
    '''
    Benjamini-Hochberg p-value correction for multiple hypothesis testing.
    adapted from https://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python/21739593

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
        np.ndarray descending_p_order, in_order, steps, adjusted_p

    pvalue = np.array(pvalue, dtype='float')
    descending_p_order = pvalue.argsort()[::-1]
    in_order = descending_p_order.argsort()
    steps = float(len(pvalue)) / np.arange(len(pvalue), 0, -1)
    adjusted_p = np.minimum(1, np.minimum.accumulate(steps * pvalue[descending_p_order]))
    return adjusted_p[in_order]