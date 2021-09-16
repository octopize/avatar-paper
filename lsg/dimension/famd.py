import sys
from itertools import chain, repeat

import numpy as np
import pandas as pd
from scipy import linalg

from lsg.dimension.dimension_reduction import DimensionReduction
from lsg.dimension.svd import SVD


class FAMD(DimensionReduction):
    """Factor Analysis of Mixed Data.

    Arguments:
        nf {int} -- number of components to keep

    Attributes:
        nf {int} -- number of components to keep
        is_quanti {array} -- list of numerical variables
        is_quali {array} -- list of categorical variables
        mean {array} -- empirical mean of each numerical variable
        std {array} -- empirical standard deviation of each numerical variable
        prop {array} -- proportion of modalities of each categorical variable
        explained_variance {array} -- amount of variance explained by the components
        explained_variance_ratio {array} -- percentage of variance explained by the components
        s {array} -- singular values corresponding to each of the selected components
        U, V {array} -- unitary matrices from the Singular Value Decomposition
        columns {array} -- columns names in the projection
    """

    def fit(self, df):
        """Fit the model with df.

        Arguments:
            df {array-like} -- training data

        Returns:
            object -- the instance itself
        """
        super().fit(df)

        if not isinstance(df, pd.DataFrame):
            df = pd.DataFrame(df)

        # select the categorical and continuous columns
        self.is_quanti = df.select_dtypes(
            include=["int", "float", "number"]
        ).columns.values
        self.is_quali = df.select_dtypes(
            exclude=["int", "float", "number"]
        ).columns.values

        # set the columns and row weights
        weight_df = pd.DataFrame([self.col_w], columns=df.columns)
        weight_quanti = weight_df[self.is_quanti]
        weight_quali = weight_df[self.is_quali]
        row_w = [1 / len(df) for i in range(len(df))]

        # scale the continuous data
        df_quanti = df[self.is_quanti]
        self.mean = np.mean(df_quanti, axis=0)
        df_quanti -= self.mean
        self.std = np.std(df_quanti, axis=0)
        self.std[self.std <= sys.float_info.min] = 1
        df_quanti /= self.std

        # get the number of modality for each quali variable
        modality_numbers = []
        for column in weight_quali.columns:
            modality_numbers += [len(df[column].unique())]

        # set weight vector for categorical columns
        weight_quali_rep = list(
            chain.from_iterable(
                repeat(i, j)
                for i, j in zip(list(weight_quali.iloc[0]), modality_numbers)
            )
        )

        # scale the categorical data
        df_quali = pd.get_dummies(df[self.is_quali].astype("category"))
        self.prop = np.mean(df_quali, axis=0)
        df_quali -= self.prop
        df_quali /= np.sqrt(self.prop)
        self.__modalities = df_quali.columns.values

        df_scale = pd.concat([df_quanti, df_quali], axis=1)

        col_w = list(weight_quanti.iloc[0]) + weight_quali_rep

        df_array = df_scale.values

        # apply the weights
        Z = ((df_array * col_w).T * row_w).T

        # compute the svd
        U, s, V = SVD(Z)
        U = ((U.T) / np.sqrt(row_w)).T
        V = V / np.sqrt(col_w)

        # compute eigenvalues and explained variance
        self.explained_variance = (s ** 2) / (df.shape[0] - 1)
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
            self.contrib, self.cos2 = self._FAMD_contrib(df, df_array, col_w, row_w)

        return self

    def fit_transform(self, df):
        """Fit the model with df and apply the dimensionality reduction on df.

        Arguments:
            df {array-like} -- training data

        Returns:
            array-like -- the transformed values
        """
        self.fit(df)

        return pd.DataFrame(self.U * self.s, columns=self.columns)

    def transform(self, df):
        """Apply dimensionality reduction to df.

        Arguments:
            df {array-like} -- data to project

        Returns:
            array-like -- the transformed values
        """
        super().transform(df)

        if not isinstance(df, pd.DataFrame):
            df = pd.DataFrame(df)

        df_quanti = df[self.is_quanti]
        df_quanti = (df_quanti - self.mean) / self.std

        # scale
        df_quali = pd.get_dummies(df[self.is_quali].astype("category"))
        for mod in self.__modalities:
            if mod not in df_quali:
                df_quali[mod] = 0
        df_quali = df_quali[self.__modalities]
        df_quali = (df_quali - self.prop) / np.sqrt(self.prop)

        df_scale = pd.concat([df_quanti, df_quali], axis=1)
        return pd.DataFrame(np.dot(df_scale, self.V.T), columns=self.columns)

    def _FAMD_contrib(self, df_original, df, col_w, row_w):
        """Compute the contribution of the variables in the axis for the FAMD.

        Arguments:
        df_orginal {dataframe} -- original dataframe
        df {ndarray} -- scaled data with dummies
        col_w {list} -- column weights
        row_w {list} -- row weights

        Returns:
            {ndarray} -- contributions of the variables
            {ndarray} -- cos2 of the variables
        """
        df2 = np.array(pd.DataFrame(df).applymap(lambda x: x ** 2))
        df = pd.DataFrame(df)

        # svd of x with row_w and col_w
        weightedTc = self._rmultiplication(
            self._rmultiplication(df.T, np.sqrt(col_w)).T, np.sqrt(row_w)
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
        dfrow_w = ((df2.T) * row_w).T
        dist2 = []
        for i in range(len(dfrow_w[0])):
            dist2 += [np.sum(dfrow_w[:, i])]
            if abs(abs(dist2[i]) - 1) < 0.001:
                dist2[i] = 1

        cor = ((coord_var.T) / np.sqrt(dist2)).T
        cos2 = cor ** 2

        # compute eta2
        df_original.index = range(len(df_original))
        dfquali = df_original[self.is_quali]
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
