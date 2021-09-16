import numpy as np
from sklearn.exceptions import NotFittedError


class DimensionReduction:
    """Base class for dimension reduction methods.

    Warning: This class should not be used directly. Use PCA, MCA or FAMD instead.

    Arguments:
        nf {int} -- number of components to keep

    Attributes:
        nf {int} -- number of components to keep
    """

    def __init__(self, nf=None, col_w=None, stats=None):
        self.nf = nf
        self.col_w = col_w
        self.stats = stats
        self.__is_fitted = False

    def fit(self, x):
        """Fit the model with x.

        Arguments:
            x {array} training data
        """
        if self.nf is None:
            self.nf = min(x.shape)
        elif self.nf <= 0:
            raise ValueError("nf", "The number of components must be positive.")

        if self.col_w is None:
            self.col_w = np.ones(x.shape[1])
        elif len(self.col_w) != x.shape[1]:
            raise ValueError(
                "col_w",
                "The weight parameter should be of size  " + str(x.shape[1]) + ".",
            )

        self.__is_fitted = True

    def fit_transform(self, x):
        """Fit the model with x and apply the dimensionality reduction on x.

        Arguments:
            x {array-like} -- training data
        """

    def transform(self, x):
        """Apply dimensionality reduction to x.

        Arguments:
            x {array-like} -- data to project
        """
        if not self.__is_fitted:
            raise NotFittedError(
                "This {} instance is not fitted yet. Call 'fit' or 'fit_transform' "
                "before using this estimator.".format(type(self).__name__)
            )
