#!/usr/bin/env python

from statsmodels.sandbox.stats.multicomp import multipletests

def padj(pvalues,**kwargs):
     '''
     This is just a re-wrap of statsmodel multipletests to make it easier to use
     '''
     rej, padj, _a , _b = multipletests(pvalues, **kwargs)
     return padj
