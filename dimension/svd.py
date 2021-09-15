from scipy import linalg
from sklearn.utils import extmath


def SVD(X):
    """Singular Value Decomposition.

    Arguments:
        X {array-like} -- matrix to decompose

    Returns:
        array -- unitary matrix having left singular vectors as columns
        array -- the singular values
        array -- unitary matrix having right singular vectors as rows
    """
    U, s, V = linalg.svd(X, full_matrices=False)
    # U, s, V = extmath.randomized_svd(X, n_components=ncp) # Truncated randomized SVD
    U, V = extmath.svd_flip(U, V)

    return U, s, V
