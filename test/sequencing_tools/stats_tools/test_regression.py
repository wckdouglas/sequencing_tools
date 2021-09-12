import numpy as np
from hypothesis import given, strategies as st
from sequencing_tools.stats_tools.regression import GradientDescent


def test_GradientDescent():
    # from sequencing_tools.stats_tools.regression import Bootstrap, GradientDescent
    # import numpy as np
    np.random.seed(123)

    size = 100
    beta = np.array([3, 4])
    X = 5 * np.random.rand(size, 2)
    y = np.matmul(X, beta) + np.random.randn(size)
    gd = GradientDescent(
        verbose=True, max_iter=100000, lr=1e-2, method="mean", limit=1e-6
    )
    gd.fit(X, y)
    assert np.allclose(gd.B, beta, atol=0.1)
