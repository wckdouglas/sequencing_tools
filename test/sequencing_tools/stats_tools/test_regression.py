import numpy as np
from hypothesis import given, strategies as st
from sequencing_tools.stats_tools.regression import GradientDescent


def test_regression():
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


@given(
    beta1=st.floats(
        min_value=0,
        max_value=20,
        exclude_min=False,
        exclude_max=False,
        allow_infinity=False,
    ),
    beta2=st.floats(
        min_value=-20,
        max_value=0,
        exclude_min=False,
        exclude_max=False,
        allow_infinity=False,
    ),
    beta3=st.floats(
        min_value=0,
        max_value=20,
        exclude_min=False,
        exclude_max=False,
        allow_infinity=False,
    ),
)
def test_regression2(beta1, beta2, beta3):
    beta = np.array([beta1, beta2, beta3])
    X = 5 * np.random.rand(100, len(beta))  # times 5 for signal-to-noise ratio
    y = np.matmul(X, beta) + np.random.randn(100)  # the last element is noise
    gd = GradientDescent(
        verbose=True, max_iter=100000, lr=1e-2, method="mean", limit=1e-6
    )
    gd.fit(X, y)
    assert np.allclose(gd.B, beta, atol=0.1)