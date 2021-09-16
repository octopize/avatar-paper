import sys

import numpy as np
import pandas as pd

from lsg.dimension.dimension_reduction import DimensionReduction
from lsg.dimension.svd import SVD


class PCA(DimensionReduction):
    """Principal Component Analysis.

    Arguments:
        nf {int} -- number of components to keep

    Attributes:
        nf {int} -- number of components to keep
        __is_fitted {bool} -- specify if the model is fitted or not
        mean {array} -- empirical mean of each variable
        std {array} -- empirical standard deviation of each variable
        explained_variance {array} -- amount of variance explained by each components
        explained_variance_ratio {array} -- percentage of variance explained by each components
        s {array} -- singular values corresponding to each of the selected components
        U, V {array} -- unitary matrices from the Singular Value Decomposition
        columns {array} -- columns names in the projection
    """

    def __init__(self, nf=None, col_w=None):
        super().__init__(nf, col_w)

    def fit(self, df, scale=True):
        """Fit the model with df.

        Arguments:
            df {array-like} -- training data
            scale {bool} -- whether or not to scale the data (default: {True})

        Returns:
            object -- the instance itself
        """
        super().fit(df)
        df = np.array(df, copy=True, dtype="float64")

        # set row weights
        row_w = [1 / len(df) for i in range(len(df))]

        # scale data
        self.mean = np.mean(df, axis=0)
        df -= self.mean
        if scale:
            self.std = np.std(df, axis=0)
            self.std[self.std <= sys.float_info.min] = 1
            df /= self.std

        # apply weights and compute svd
        Z = ((df * self.col_w).T * row_w).T
        U, s, V = SVD(Z)

        U = ((U.T) / np.sqrt(row_w)).T
        V = V / np.sqrt(self.col_w)

        # compute eigenvalues and explained variance
        self.explained_variance = (s ** 2) / (df.shape[0] - 1)
        self.explained_variance_ratio = (
            self.explained_variance / self.explained_variance.sum()
        )[: self.nf]
        self.explained_variance = self.explained_variance[: self.nf]

        self.U = U[:, : self.nf]
        self.s = s[: self.nf]
        self.V = V[: self.nf, :]

        self.columns = ["Dim. {}".format(i + 1) for i in range(self.nf)]
        return self

    def fit_transform(self, df, scale=True):
        """Fit the model with df and apply the dimensionality reduction on df.

        Arguments:
            df {array-like} -- training data
            scale {bool} -- whether or not to scale the data (default: {True})

        Returns:
            array-like -- the transformed values
        """
        self.fit(df, scale=scale)

        return pd.DataFrame(self.U * self.s, columns=self.columns)

    def transform(self, df):
        """Apply dimensionality reduction to df.

        Arguments:
            df {array-like} -- data to project

        Returns:
            array-like -- the transformed values
        """
        super().transform(df)
        df = np.array(df, copy=True, dtype="float64")

        # scale
        df -= self.mean
        df /= self.std
        return pd.DataFrame(np.dot(df, self.V.T), columns=self.columns)
