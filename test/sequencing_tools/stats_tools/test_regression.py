import numpy as np
from sequencing_tools.stats_tools.regression import GradientDescent


def test_regression():
    # from sequencing_tools.stats_tools.regression import Bootstrap, GradientDescent
    # import numpy as np
    np.random.seed(123)

    X = 5 * np.random.rand(100, 2)
    y = np.matmul(X, np.array([3, 4])) + np.random.randn(100)
    gd = GradientDescent(
        verbose=True, max_iter=100000, lr=1e-2, method="mean", limit=1e-6
    )
    gd.fit(X, y)
    assert np.abs(gd.B[0] - 3) < 0.1
    assert np.abs(gd.B[1] - 4) < 0.1
