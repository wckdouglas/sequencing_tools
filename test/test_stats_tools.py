from sequencing_tools.stats_tools import p_adjust, levenshtein_distance, hamming_distance
from sequencing_tools.stats_tools.regression import Bootstrap, GradientDescent
import numpy as np

def test_padjust():
    '''
    in R:
    test_p <- c(  2.49623261e-03,   2.43927089e-02,   8.80158825e-03,
         1.09798521e-05,   8.14573043e-01,   4.20958895e-01,
         3.23475673e-01,   1.65218466e-01,   3.87039350e-02,
         4.81924429e-03,   7.15204680e-04,   2.37862808e-01,
         2.77650933e-02,   2.25735871e-02,   6.96230468e-03,
         3.13598359e-03,   2.84268503e-03,   1.46484891e-03,
         1.20357262e-03,   7.83426751e-01,   2.23548765e-01,
         4.12654807e-02,   2.72907764e-02,   2.55000319e-02,
         1.82854508e-02,   1.45822194e-02,   1.34166230e-02,
         8.39874556e-01,   5.28877785e-01,   4.85792040e-01,
         4.52649675e-01,   3.50140056e-01,   3.27987681e-01,
         2.14560873e-01,   2.07004665e-01,   2.03179557e-01,
         1.87989254e-01,   1.73866175e-01,   1.72914799e-01,
         1.61658121e-01,   1.57143261e-01,   1.55072731e-01,
         1.52475888e-01,   1.45224636e-01,   1.39084689e-01,
         1.22256279e-01,   1.20853414e-01,   1.11466566e-01,
         1.09664995e-01,   1.05955886e-01)
    round(p.adjust(test_p, method='BH'),3)
    '''

    test_p = np.array([2.354054**-7,2.101590**-5,2.576842**-5,9.814783**-5,1.052610**-4,1.241481**-4,
        1.325988**-4,1.568503**-4,2.254557**-4,3.795380**-4,6.114943**-4,1.613954**-3,
        3.302430**-3,3.538342**-3,5.236997**-3,6.831909**-3,7.059226**-3,8.805129**-3,
        9.401040**-3,1.129798**-2,2.115017**-2,4.922736**-2,6.053298**-2,6.262239**-2,
        7.395153**-2,8.281103**-2,8.633331**-2,1.190654**-1,1.890796**-1,2.058494**-1,
        2.209214**-1,2.856000**-1,3.048895**-1,4.660682**-1,4.830809**-1,4.921755**-1,
        5.319453**-1,5.751550**-1,5.783195**-1,6.185894**-1,6.363620**-1,6.448587**-1,
        6.558414**-1,6.885884**-1,7.189864**-1,8.179539**-1,8.274487**-1,8.971300**-1,
        9.118680**-1,9.437890**-1])

    expected_p = np.array([ 0.022,  0.077,  0.044,  0.001,  0.831,  0.478,  0.39 ,  0.256,
        0.102,  0.03 ,  0.018,  0.297,  0.077,  0.077,  0.039,  0.022,
        0.022,  0.018,  0.018,  0.816,  0.287,  0.103,  0.077,  0.077,
        0.07 ,  0.061,  0.061,  0.84 ,  0.563,  0.528,  0.503,  0.407,
        0.39 ,  0.282,  0.28 ,  0.28 ,  0.269,  0.256,  0.256,  0.256,
        0.256,  0.256,  0.256,  0.256,  0.256,  0.245,  0.245,  0.242,
        0.242,  0.242])

    padj = p_adjust(test_p)
    assert(np.allclose(np.round(padj,3), expected_p))


def test_distance():
    assert(levenshtein_distance('AACCA','AACCT') == 1)
    assert(hamming_distance('AACCA','AACCT') == 1)


def test_bootstrap():
    n = 100
    group_size = 10 
    a = np.random.rand(100)
    bs = Bootstrap()
    for i, b in enumerate(bs.bootstrap(a, group_size=group_size, n_boots=n)):
        assert( len(b) == group_size)
    
    assert( i == n-1 )
    

def test_regression():
    #from sequencing_tools.stats_tools.regression import Bootstrap, GradientDescent
    #import numpy as np
    np.random.seed(123)
    X = 2 * np.random.rand(100,2) 
    y = np.sum(np.array([3,4]) * X, axis=1)+np.random.randn(100) 

    gd = GradientDescent(verbose = True, max_iter=100000, lr = 1e-2, method='mean')
    gd.fit(X,y)
    assert np.abs(gd.B[0] - 3.5) < 0.2
    assert np.abs(gd.B[1] - 3.2) < 0.2