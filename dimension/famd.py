import numpy as np
import pandas as pd
import sys
from itertools import repeat, chain
import numpy as np
from scipy import linalg

from dimension_reduction import DimensionReduction
from svd import SVD


class FAMD(DimensionReduction):
    """Factor Analysis of Mixed Data

    Arguments:
        nf {int} -- number of components to keep

    Attributes:
        nf {int} -- number of components to keep
        __is_fitted {bool} -- specify if the model is fitted or not
        is_quanti {array} -- list of numerical variables
        is_quali {array} -- list of categorical variables
        mean {array} -- empirical mean of each numerical variable
        std {array} -- empirical standard deviation of each numerical variable
        prop {array} -- proportion of modalities of each categorical variable
        __modalities {array} -- list of modalities of each categorical variables
        explained_variance {array} -- amount of variance explained by each of the selected components
        explained_variance_ratio {array} -- percentage of variance explained by each of the selected components
        s {array} -- singular values corresponding to each of the selected components
        U, V {array} -- unitary matrices from the Singular Value Decomposition
        columns {array} -- columns names in the projection
    """

    def __init__(self, nf=None, col_w=None, stats=True):
        super().__init__(nf, col_w, stats)

    def fit(self, X):
        """Fit the model with X.

        Arguments:
            X {array-like} -- training data

        Returns:
            object -- the instance itself
        """
        super().fit(X)

        if not isinstance(X, pd.DataFrame):
            X = pd.DataFrame(X)

        # select the categorical and continuous columns
        self.is_quanti = X.select_dtypes(
            include=["int", "float", "number"]
        ).columns.values
        self.is_quali = X.select_dtypes(
            exclude=["int", "float", "number"]
        ).columns.values

        # set the columns and row weights
        weight_df = pd.DataFrame([self.col_w], columns=X.columns)
        weight_quanti = weight_df[self.is_quanti]
        weight_quali = weight_df[self.is_quali]
        row_w = [1 / len(X) for i in range(len(X))]

        # scale the continuous data
        X_quanti = X[self.is_quanti]
        self.mean = np.mean(X_quanti, axis=0)
        X_quanti -= self.mean
        self.std = np.std(X_quanti, axis=0)
        self.std[self.std <= sys.float_info.min] = 1
        X_quanti /= self.std

        # get the number of modality for each quali variable
        modality_numbers = []
        for column in weight_quali.columns:
            modality_numbers += [len(X[column].unique())]

        # set weight vector for categorical columns
        weight_quali_rep = list(
            chain.from_iterable(
                repeat(i, j)
                for i, j in zip(list(weight_quali.iloc[0]), modality_numbers)
            )
        )

        # scale the categorical data
        X_quali = pd.get_dummies(X[self.is_quali].astype("category"))
        self.prop = np.mean(X_quali, axis=0)
        X_quali -= self.prop
        X_quali /= np.sqrt(self.prop)
        self.__modalities = X_quali.columns.values

        X_scale = pd.concat([X_quanti, X_quali], axis=1)

        col_w = list(weight_quanti.iloc[0]) + weight_quali_rep

        X_array = X_scale.values

        # apply the weights
        Z = ((X_array * col_w).T * row_w).T

        # compute the svd
        U, s, V = SVD(Z)
        U = ((U.T) / np.sqrt(row_w)).T
        V = V / np.sqrt(col_w)

        # compute eigenvalues and explained variance
        self.explained_variance = (s ** 2) / (X.shape[0] - 1)
        self.explained_variance_ratio = (
            self.explained_variance / self.explained_variance.sum()
        )[: self.nf]
        self.explained_variance = self.explained_variance[: self.nf]

        self.U = U[:, : self.nf]
        self.s = s[: self.nf]
        self.V = V[: self.nf, :]

        # compute contributions and cos2
        self.columns = ["Dim. {}".format(i + 1) for i in range(len(self.s))]
        if self.stats:
            self.contrib, self.cos2 = self._FAMD_contrib(X, X_array, col_w, row_w)

        return self

    def fit_transform(self, X):
        """Fit the model with X and apply the dimensionality reduction on X.

        Arguments:
            X {array-like} -- training data

        Returns:
            array-like -- the transformed values
        """
        self.fit(X)

        return pd.DataFrame(self.U * self.s, columns=self.columns)

    def transform(self, X):
        """Apply dimensionality reduction to X.

        Arguments:
            X {array-like} -- data to project

        Returns:
            array-like -- the transformed values
        """
        super().transform(X)

        if not isinstance(X, pd.DataFrame):
            X = pd.DataFrame(X)

        X_quanti = X[self.is_quanti]
        X_quanti = (X_quanti - self.mean) / self.std

        # scale
        X_quali = pd.get_dummies(X[self.is_quali].astype("category"))
        for mod in self.__modalities:
            if mod not in X_quali:
                X_quali[mod] = 0
        X_quali = X_quali[self.__modalities]
        X_quali = (X_quali - self.prop) / np.sqrt(self.prop)

        X_scale = pd.concat([X_quanti, X_quali], axis=1)
        return pd.DataFrame(np.dot(X_scale, self.V.T), columns=self.columns)

    def _FAMD_contrib(self, X_original, X, col_w, row_w):
        """Compute the contribution of the variables in the axis for the FAMD

        Arguments :
        X_orginal {dataframe} -- original dataframe
        X {ndarray} -- scaled data with dummies
        col_w {list} -- column weights
        row_w {list} -- row weights

        Returns :
            {ndarray} -- contributions of the variables
            {ndarray} -- cos2 of the variables
        """
        X2 = np.array(pd.DataFrame(X).applymap(lambda x: x ** 2))
        X = pd.DataFrame(X)

        # svd of x with row_w and col_w
        weightedTc = self._rmultiplication(
            self._rmultiplication(X.T, np.sqrt(col_w)).T, np.sqrt(row_w)
        )
        U, s, V = linalg.svd(weightedTc.T, full_matrices=False)
        ncp0 = min(len(weightedTc.iloc[0]), len(weightedTc), self.nf)
        U = U[:, :ncp0]
        V = V.T[:, :ncp0]
        s = s[:ncp0]
        tmp = V
        V = U
        U = tmp
        mult = np.sign(np.sum(V, axis=0))

        # final V
        mult1 = pd.DataFrame(
            np.array(
                pd.DataFrame(np.array(self._rmultiplication(pd.DataFrame(V.T), mult)))
            ).T
        )
        V = pd.DataFrame()
        for i in range(len(mult1)):
            V[i] = mult1.iloc[i] / np.sqrt(col_w[i])
        V = np.array(V).T

        # final U
        mult1 = pd.DataFrame(
            np.array(
                pd.DataFrame(np.array(self._rmultiplication(pd.DataFrame(U.T), mult)))
            ).T
        )
        U = pd.DataFrame()
        for i in range(len(mult1)):
            U[i] = mult1.iloc[i] / np.sqrt(row_w[i])
        U = np.array(U).T
        eig = s ** 2
        # end of the svd

        # compute the contribution
        coord_var = np.array(V[0] * s)
        for i in range(1, len(V[:, 0])):
            coord_var = np.vstack((coord_var, V[i] * s))
        contrib_var = (((((coord_var ** 2) / eig).T) * col_w).T) * 100

        # compute cos2
        Xrow_w = ((X2.T) * row_w).T
        dist2 = []
        for i in range(len(Xrow_w[0])):
            dist2 += [np.sum(Xrow_w[:, i])]
            if abs(abs(dist2[i]) - 1) < 0.001:
                dist2[i] = 1

        cor = ((coord_var.T) / np.sqrt(dist2)).T
        cos2 = cor ** 2

        # compute eta2
        X_original.index = range(len(X_original))
        dfquali = X_original[self.is_quali]
        eta2 = []
        fi = 0
        coord = pd.DataFrame(
            self.U[:, :ncp0] * self.s[:ncp0], columns=self.columns[:ncp0]
        )
        mods = []
        # for each qualitative column in the original data set
        for count, col in enumerate(dfquali.columns):
            dummy = pd.get_dummies(dfquali[col].astype("category"))
            mods += [len(dummy.columns) - 1]
            # for each dimension
            dim = []
            for j, coordcol in enumerate(coord.columns):
                # for each modality of the qualitative column
                s = 0
                for i in range(len(dummy.columns)):
                    s += (
                        np.array(dummy.T)[i] * coord[coordcol] * row_w
                    ).sum() ** 2 / self.prop[fi + i]
                dim += [s]
            eta1 = (
                np.array(dim) / (np.array((coord ** 2)).T * row_w).sum(axis=1)
            ).tolist()
            eta2 += [eta1]
            fi += len(dummy.columns)

            cos2 = cos2[: len(self.is_quanti)]

        cos2 = cos2 ** 2
        eta2 = np.array(eta2) ** 2
        eta2 = (eta2.T / mods).T

        cos2 = np.concatenate([cos2, eta2], axis=0)

        return contrib_var, cos2

    @staticmethod
    def _rmultiplication(F, marge):
        multiplication = pd.DataFrame()
        for col in F.columns:
            multiplication[col] = F[col] * marge
        multiplication.index = F.index
        return multiplication
