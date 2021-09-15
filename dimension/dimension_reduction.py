from sklearn.exceptions import NotFittedError
import numpy as np


class DimensionReduction:
    """Base class for dimension reduction methods.

    Warning: This class should not be used directly. Use PCA, MCA or FAMD instead.

    Arguments:
        nf {int} -- number of components to keep

    Attributes:
        nf {int} -- number of components to keep
        __is_fitted {bool} -- specify if the model is fitted or not
    """

    def __init__(self, nf=None, col_w=None, stats=None):
        self.nf = nf
        self.col_w = col_w
        self.stats = stats
        self.__is_fitted = False

    def fit(self, X):
        """Fit the model with X.

        Arguments:
            X {array-like} -- training data
        """
        if self.nf is None:
            self.nf = min(X.shape)
        elif self.nf <= 0:
            raise ValueError("nf", "The number of components must be positive.")

        if self.col_w is None:
            self.col_w = np.ones(X.shape[1])
        elif len(self.col_w) != X.shape[1]:
            raise ValueError(
                "col_w",
                "The weight parameter should be of size  " + str(X.shape[1]) + ".",
            )

        self.__is_fitted = True

    def fit_transform(self, X):
        """Fit the model with X and apply the dimensionality reduction on X.

        Arguments:
            X {array-like} -- training data
        """

    def transform(self, X):
        """Apply dimensionality reduction to X.

        Arguments:
            X {array-like} -- data to project
        """
        if not self.__is_fitted:
            raise NotFittedError(
                "This {} instance is not fitted yet. Call 'fit' or 'fit_transform' before using this estimator.".format(
                    type(self).__name__
                )
            )
