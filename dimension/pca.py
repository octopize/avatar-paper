import numpy as np
import pandas as pd
import sys

from dimension.dimension_reduction import DimensionReduction
from dimension.svd import SVD


class PCA(DimensionReduction):
    """Principal Component Analysis

    Arguments:
        nf {int} -- number of components to keep

    Attributes:
        nf {int} -- number of components to keep
        __is_fitted {bool} -- specify if the model is fitted or not
        mean {array} -- empirical mean of each variable
        std {array} -- empirical standard deviation of each variable
        explained_variance {array} -- amount of variance explained by each of the selected components
        explained_variance_ratio {array} -- percentage of variance explained by each of the selected components
        s {array} -- singular values corresponding to each of the selected components
        U, V {array} -- unitary matrices from the Singular Value Decomposition
        columns {array} -- columns names in the projection
    """

    def __init__(self, nf=None, col_w=None):
        super().__init__(nf, col_w)

    def fit(self, X, scale=True):
        """Fit the model with X.

        Arguments:
            X {array-like} -- training data
            scale {bool} -- whether or not to scale the data (default: {True})

        Returns:
            object -- the instance itself
        """
        super().fit(X)
        X = np.array(X, copy=True, dtype="float64")

        # set row weights
        row_w = [1 / len(X) for i in range(len(X))]

        # scale data
        self.mean = np.mean(X, axis=0)
        X -= self.mean
        if scale:
            self.std = np.std(X, axis=0)
            self.std[self.std <= sys.float_info.min] = 1
            X /= self.std

        # apply weights and compute svd
        Z = ((X * self.col_w).T * row_w).T
        U, s, V = SVD(Z)

        U = ((U.T) / np.sqrt(row_w)).T
        V = V / np.sqrt(self.col_w)

        # compute eigenvalues and explained variance
        self.explained_variance = (s ** 2) / (X.shape[0] - 1)
        self.explained_variance_ratio = (
            self.explained_variance / self.explained_variance.sum()
        )[: self.nf]
        self.explained_variance = self.explained_variance[: self.nf]

        self.U = U[:, : self.nf]
        self.s = s[: self.nf]
        self.V = V[: self.nf, :]

        self.columns = ["Dim. {}".format(i + 1) for i in range(self.nf)]
        return self

    def fit_transform(self, X, scale=True):
        """Fit the model with X and apply the dimensionality reduction on X.

        Arguments:
            X {array-like} -- training data
            scale {bool} -- whether or not to scale the data (default: {True})

        Returns:
            array-like -- the transformed values
        """
        self.fit(X, scale=scale)

        return pd.DataFrame(self.U * self.s, columns=self.columns)

    def transform(self, X):
        """Apply dimensionality reduction to X.

        Arguments:
            X {array-like} -- data to project

        Returns:
            array-like -- the transformed values
        """
        super().transform(X)
        X = np.array(X, copy=True, dtype="float64")

        # scale
        X -= self.mean
        X /= self.std
        return pd.DataFrame(np.dot(X, self.V.T), columns=self.columns)
