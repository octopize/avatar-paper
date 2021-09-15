from dimension.pca import PCA
from dimension.mca import MCA
from dimension.famd import FAMD
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


class Projection:
    def __init__(self):
        pass

    def _variable_correlation(self):
        """Compute the correlation between the axis' and the variables

        Returns :
            Correlation matrix
        """

        # select columns and project data
        self._df_quanti = self._data_frame[self._is_quanti]
        coord = self.model.transform(self._data_frame)

        if len(self._is_quali) > 0:
            self._df_quali = pd.get_dummies(
                self._data_frame[self._is_quali].astype("category")
            )
            bind = pd.concat([self._df_quanti, self._df_quali], axis=1)
        else:
            bind = self._df_quanti

        cor = pd.DataFrame(
            {
                component: {
                    feature: coord[component].corr(bind[feature].reset_index(drop=True))
                    for feature in bind.columns
                }
                for component in coord.columns
            }
        )
        return cor

    def fit_transform(self, data_frame, nf=None, col_w=None, stats=True):
        """Project data into a lower dimensional space using PCA,MCA or FAMD.

        Arguments:
            data_frame {array-like} -- matrix to project
            nf {int} -- number of components to keep (default: {min(data_frame.shape[0], data_frame.shape[1])})
            col_w {numeric list} -- allow to increase the importance of a variable to separate records on the first PC's

        Returns:
            array-like -- the transformed values
            object -- the dimension reduction model
        """
        # Transform date into time post epoch
        datetime_variables = []
        for i in range(0, data_frame.shape[1]):
            if data_frame.iloc[:, i].dtype == ("datetime64[ns]"):
                data_frame.iloc[:, i] = (
                    data_frame.iloc[:, i] - np.datetime64("1970-01-01T00:00:00Z")
                ) / np.timedelta64(1, "s")
                datetime_variables.append(i)
        self.datetime_variables = datetime_variables

        # check column types
        self._is_quanti = data_frame.select_dtypes(
            include=["int", "float", "number"]
        ).columns.values
        self._is_quali = data_frame.select_dtypes(
            exclude=["int", "float", "number"]
        ).columns.values

        if nf == "all":
            nf = len(pd.get_dummies(data_frame).columns.values)
        if nf == None:
            nf = 10

        # specify the correct model
        if self._is_quali.size == 0:
            print(
                "Regarding variables types, a PCA is performed for dimension reduction"
            )
            model = PCA(nf=nf, col_w=col_w)
        elif self._is_quanti.size == 0:
            print(
                "Regarding variables types, a MCA is performed for dimension reduction"
            )
            model = MCA(nf=nf, col_w=col_w, stats=stats)
        else:
            print(
                "Regarding variables types, a FAMD is performed for dimension reduction"
            )
            model = FAMD(nf=nf, col_w=col_w, stats=stats)

        model.fit_transform(data_frame)

        self.model = model
        self._data_frame = data_frame
        self.nf = model.nf

        # compute axes coordinates
        if self._is_quanti.size == 0:
            self.variable_coord = pd.DataFrame(np.dot(model.D_c, model.V.T))
        else:
            self.variable_coord = pd.DataFrame(self.model.V.T)

        # compute projection stats
        if stats:
            self.cor = self._variable_correlation()
            self.variable_coord.columns = self.cor.columns
            self.variable_coord.index = list(self.cor.index)

            if self._is_quali.size == 0:
                self.cos2 = self.cor.applymap(lambda x: x ** 2)
                self.contrib = self.cos2.div(self.cos2.sum(axis=0), axis=1).applymap(
                    lambda x: x * 100
                )
            elif self._is_quanti.size == 0:
                self.cos2 = self.cor.applymap(lambda x: x ** 2)
                self.contrib = pd.DataFrame(
                    self.model.contrib,
                    columns=self.cor.columns,
                    index=list(self.cor.index),
                )
            else:
                self.cos2 = pd.DataFrame(
                    self.model.cos2, index=list(model.is_quanti) + list(model.is_quali)
                )
                self.contrib = pd.DataFrame(
                    self.model.contrib,
                    columns=self.cor.columns,
                    index=list(self.cor.index),
                )

        return model.transform(data_frame), model

    def transform(self, X):
        """Projects the data into the fitted space

        Arguments :
            X {dataframe} -- data to project

        Raises :
            RuntimeError: 'Class should be fitted before plotting'

        Returns :
            {dataframe} -- projected data
        """

        # Transform date into time post epoch
        for i in self.datetime_variables:
            X.iloc[:, i] = (
                X.iloc[:, i] - np.datetime64("1970-01-01T00:00:00Z")
            ) / np.timedelta64(1, "s")

        if not hasattr(self, "model"):
            raise RuntimeError("Class should be fitted before plotting")
        return self.model.transform(X)

    def inverse_transform(self, coord, shuffle=False):
        """Compute the inverse projection of the coordinates

        Arguments :
            coord {dataframe} -- coordinate to inverse
            shuffle {boolean} -- False to keep similar decimals as original data set

        Returns :
            {dataframe} -- inversed data
        """
        model = self.model

        # if PCA or FAMD conpute the continuous variables
        if len(self._is_quanti) != 0:

            X = np.dot(coord, self.variable_coord.T)
            X_quanti = X[:, : len(self._is_quanti)]

            # descale
            std = np.array(model.std)
            mean = np.array(model.mean)
            inverse_quanti = (X_quanti * std) + mean
            inverse_quanti = pd.DataFrame(inverse_quanti, columns=list(self._is_quanti))

            # round to the right decimal
            for column in inverse_quanti.columns:
                inverse_quanti["decimals"] = self._data_frame[column].apply(
                    self.decimal_count
                )
                # shuffling the decimals for the avatarization
                if shuffle:
                    inverse_quanti["decimals"] = np.random.permutation(
                        inverse_quanti["decimals"].values
                    )

                inverse_quanti[column] = inverse_quanti[[column, "decimals"]].apply(
                    lambda x: np.round(x[column], int(x["decimals"])), axis=1
                )
                inverse_quanti.drop(["decimals"], axis=1, inplace=True)

            # if FAMD descale the categorical variables
            if len(self._is_quali) != 0:
                X_quali = X[:, len(self._is_quanti) :]
                prop = np.array(model.prop)
                X_quali = (X_quali) * (np.sqrt(prop)) + prop

        # if MCA no descaling
        else:
            X_quali = np.dot(coord, np.dot(model.D_c, model.V.T).T)

        # compute the categorical variables
        if len(self._is_quali) != 0:
            inverse_quali = pd.DataFrame()
            X_quali = pd.DataFrame(X_quali)
            X_quali.columns = list(
                pd.get_dummies(
                    self._data_frame[self._is_quali],
                    prefix=["" for i in range(len(self._is_quali))],
                    prefix_sep="",
                ).columns
            )

            modalities = []
            for column in self._data_frame[self._is_quali].columns:
                modalities += [len(self._data_frame[column].unique())]
            val = 0
            for i in range(len(modalities)):
                inverse_quali[
                    list(self._data_frame[self._is_quali].columns)[i]
                ] = X_quali.iloc[:, val : val + modalities[i]].idxmax(1)
                val += modalities[i]
            inverse_quali = inverse_quali.astype("str")

        # concatenate the continuous and categorical
        if len(self._is_quali) != 0 and len(self._is_quanti) != 0:
            inverse = pd.concat([inverse_quali, inverse_quanti], axis=1)
        elif len(self._is_quanti) != 0:
            inverse = inverse_quanti
        else:
            inverse = inverse_quali

        # Cast columns to same type as input
        for column in self._data_frame.columns:
            column_type = self._data_frame.loc[:, column].dtype
            inverse[column] = inverse[column].astype(column_type)

        # reorder back columns
        inverse = inverse[self._data_frame.columns]

        # Turn back datetime variables to original dtype
        for i in self.datetime_variables:
            inverse.iloc[:, i] = (
                inverse.iloc[:, i] * np.timedelta64(1, "s")
            ) + np.datetime64("1970-01-01T00:00:00Z")

        return inverse

    @staticmethod
    def decimal_count(number):
        f = str(number)
        if "." in f:
            digits = f[::-1].find(".")
        else:
            digits = 0
        return digits

    def plot_circle(self, dimensions=[1, 2], min_cor=0.1, max_var=7):
        """Plot correlation graph

        Arguments :
        dimensions {list} -- list of the dimensions to help by each axis

        Raises :
            RuntimeError: 'Class should be fitted before plotting'

        Returns :
            graph -- plot of the correlation circle
        """
        # Dimensions start from 1

        if not hasattr(self, "model"):
            raise RuntimeError("Class should be fitted before plotting")

        # Plotting circle
        figure_axis_size = 6
        explained_var_ratio = self.model.explained_variance_ratio

        fig_res = plt.figure(figsize=(figure_axis_size, figure_axis_size))
        circle1 = plt.Circle((0, 0), radius=1, color="k", fill=False)
        fig = plt.gcf()
        fig.gca().add_artist(circle1)

        # Order dataframe
        cor = self.cor.copy()
        cor["sum"] = cor.apply(
            lambda x: abs(x[dimensions[0] - 1]) + abs(x[dimensions[1] - 1]), axis=1
        )
        cor.sort_values(by="sum", ascending=False, inplace=True)

        # Plotting arrows
        texts = []
        i = 0
        for name, row in cor.iterrows():
            if i < max_var and (
                np.abs(row[dimensions[0] - 1]) > min_cor
                or np.abs(row[dimensions[1] - 1]) > min_cor
            ):
                x = row[dimensions[0] - 1]
                y = row[dimensions[1] - 1]
                plt.arrow(
                    0.0,
                    0.0,
                    x,
                    y,
                    color="k",
                    length_includes_head=True,
                    head_width=0.05,
                )

                plt.plot([0.0, x], [0.0, y], "k-")
                texts.append(plt.text(x, y, name, fontsize=2 * figure_axis_size))
                i += 1

        # Plotting vertical lines
        plt.plot([-1.1, 1.1], [0, 0], "k--")
        plt.plot([0, 0], [-1.1, 1.1], "k--")

        # Setting limits and title
        plt.xlim((-1.1, 1.1))
        plt.ylim((-1.1, 1.1))
        plt.title("Correlation Circle", fontsize=figure_axis_size * 3)

        plt.xlabel(
            "Dim "
            + str(dimensions[0])
            + " (%s%%)" % str(explained_var_ratio[dimensions[0] - 1] * 100)[:4],
            fontsize=figure_axis_size * 2,
        )
        plt.ylabel(
            "Dim "
            + str(dimensions[1])
            + " (%s%%)" % str(explained_var_ratio[dimensions[1] - 1] * 100)[:4],
            fontsize=figure_axis_size * 2,
        )

    def plot_var_contribution(self, dim=1, max_var=10, min_contrib=0.1):
        """Plot the variable contribution for a given dimension

        Arguments :
            dim {int} -- value of the dimension to plot

        Keyword Arguments :
            max_var {int} -- maximum number of variables to plot
            min_contrib {int} -- lower threshold for the variable contributions

        Raises :
            RuntimeError: 'Class should be fitted before plotting'

        Returns :
            graph of the contribution percentages per variables

        """
        # Dimensions start from 1

        if not hasattr(self, "contrib"):
            raise RuntimeError("Class should be fitted before plotting")

        # get the useful contributions
        var_contrib = self.contrib[self.contrib.columns[dim - 1]]
        if len(var_contrib) > max_var:
            var_contrib = var_contrib[:max_var]

        # check threshold
        var_contrib = [var for var in var_contrib if var > min_contrib]
        var_contrib = pd.DataFrame(var_contrib)[0]

        indices = list((-var_contrib).argsort())
        names = [list(self.contrib.index)[indices[i]] for i in range(len(indices))]

        # plot
        plt.figure(figsize=(12, 6))
        plt.bar(range(len(var_contrib)), var_contrib[indices], align="center")
        plt.xticks(range(len(var_contrib)), names, rotation="horizontal")

        # setting labels and title
        plt.title("Variables contributions to Dim. " + str(dim))
        plt.ylabel("Importance")
        plt.xlabel("Variables")
        plt.show()

    def plot_explained_var(self, max_dims=10):
        """Plot explained variance per dimension

        Arguments:
            max_dims {int} -- maximum number of dimensions to plot

        Raises:
            RuntimeError: 'Class should be fitted before plotting'

        Return:
            graph -- plot of the explained variance
        """

        if not hasattr(self, "model"):
            raise RuntimeError("Class should be fitted before plotting")

        explained_percentage = self.model.explained_variance_ratio * 100
        if len(explained_percentage) > max_dims:
            explained_percentage = explained_percentage[:max_dims]

        # plot
        plt.figure(figsize=(12, 6))
        plt.bar(range(len(explained_percentage)), explained_percentage, align="center")
        plt.xticks(
            range(len(explained_percentage)),
            range(1, len(explained_percentage) + 1),
            rotation="horizontal",
        )

        # setting labels and title
        plt.title("Explained variance plot")
        plt.ylabel("Percentage of explained variance")
        plt.xlabel("Dimensions")
        plt.show()
