import numpy as np
import pandas as pd
from scipy.sparse import diags
from scipy import linalg
import warnings
warnings.filterwarnings("ignore")
from itertools import repeat, chain

from dimension_reduction import DimensionReduction
from svd import SVD

class MCA(DimensionReduction):
    """ Multiple Correspondence Analysis 
    
    Arguments:
        nf {int} -- number of components to keep
    
    Attributes:
        nf {int} -- number of components to keep
        __is_fitted {bool} -- specify if the model is fitted or not
        __modalities {array} -- list of modalities of each categorical variables
        D_r, D_c {array} -- weight matrices
        explained_variance {array} -- amount of variance explained by each of the selected components
        explained_variance_ratio {array} -- percentage of variance explained by each of the selected components
        s {array} -- singular values corresponding to each of the selected components
        U, V {array} -- unitary matrices from the Singular Value Decomposition
        columns {array} -- columns names in the projection
    """
    def __init__(self, nf=None, col_w = None, stats = True):
        super().__init__(nf, col_w, stats)
    
    def __fit(self, X):
        """ Fit the model with X.

        Arguments:
            X {array-like} -- training data
        
        Returns:
            object -- the instance itself
            array-like -- the modified training data
        """
        super().fit(X)
        if not isinstance(X, pd.DataFrame):
            X = pd.DataFrame(X)
        
        X_original = X.copy()
            
        #initiate row and columns weights
        row_w = [1/len(X)for i in range(len(X))]
        modality_numbers = []
        for column in X.columns:
            modality_numbers += [len(X[column].unique())]    
        col_w = list(chain.from_iterable(repeat(i, j) for i, j in zip(self.col_w, modality_numbers)))    

        X = pd.get_dummies(X.astype("category"))
        self.__modalities = X.columns.values
        
        # scale data
        X /= X.sum().sum()
        c = np.sum(X, axis=0)
        r = np.sum(X, axis=1)
    
        # set D_r and D_c
        eps = np.finfo(float).eps
        if X.shape[0] >= 10000:
            self.D_r = diags(1/(eps+np.sqrt(r)))
        else:
            self.D_r = np.diag(1/(eps+np.sqrt(r)))
        self.D_c = np.diag(1/(eps+np.sqrt(c)))

        T = self.D_r @ (X - np.outer(r, c)) @ self.D_c
        X /= r[:, None]
        
        # apply the weights and compute the svd
        Z = ((T*col_w).T*row_w).T
        U, s, V = SVD(Z)
        
        # compute eigenvalues and explained variance
        self.explained_variance = (s**2) / (X.shape[0] - 1)
        self.explained_variance_ratio = (self.explained_variance / self.explained_variance.sum())[:self.nf]
        self.explained_variance = self.explained_variance[:self.nf]
        
        #apply similar products to U and V
        U = ((U.T)/np.sqrt(row_w)).T
        V = V/np.sqrt(col_w)
        
        self.U = U[:, :self.nf]
        self.s = s[:self.nf]
        self.V = V[:self.nf, :]

        self.columns = ["Dim. {}".format(i+1) for i in range(min(self.nf, self.s.shape[0]))]

        # compute contributions
        if self.stats:
            self.contrib = self._MCA_contrib(X_original)

        return self, X
    
    def fit(self, X):
        """ Fit the model with X.

        Arguments:
            X {array-like} -- training data
        
        Returns:
            object -- the instance itself
        """
        self.__fit(X)
        return self

    def fit_transform(self, X):
        """ Fit the model with X and apply the dimensionality reduction on X.

        Arguments:
            X {array-like} -- training data
        
        Returns:
            array-like -- the transformed values
        """
        _, X = self.__fit(X)
        return pd.DataFrame(np.dot(X, np.dot(self.D_c, self.V.T)), columns=self.columns)
    
    def transform(self, X):
        """ Apply dimensionality reduction to X.

        Arguments:
            X {array-like} -- data to project

        Returns:
            array-like -- the transformed values
        """
        super().transform(X)
        
        if not isinstance(X, pd.DataFrame):
            X = pd.DataFrame(X)
        X = pd.get_dummies(X.astype("category"))
        for mod in self.__modalities:
            if mod not in X:
                X[mod] = 0
        X = X[self.__modalities]

        # scale
        X /= X.sum().sum()
        X /= np.sum(X, axis=1)[:, None]
        
        return pd.DataFrame(np.dot(X, np.dot(self.D_c, self.V.T)), columns=self.columns)

    @staticmethod
    def _rmultiplication(F,marge):
        multiplication = pd.DataFrame()
        for col in F.columns:
            multiplication[col] = F[col]*marge
        multiplication.index = F.index
        return multiplication

    @staticmethod
    def _rdivision(F,marge):
        division = pd.DataFrame()
        for col in F.columns:
            division[col] = F[col]/marge
        division.index = F.index
        return division

    def _MCA_contrib(self, X_original):
        '''Compute the contribution of the variables in the axes for the PCA

        Arguments :
            model {object} -- MCA model

        Returns:
            {dataframe} -- contributions of the variables
        '''

        V = np.dot(self.D_c, self.V.T)
        total = pd.get_dummies(X_original.astype("category")).sum().sum()
        X = pd.get_dummies(X_original.astype("category"))
        F = X/total

        # Column and row weights
        marge_col = F.sum(axis=0)
        marge_row = F.sum(axis=1)
        fsurmargerow = self._rdivision(F,marge_row)
        fmargerowT = pd.DataFrame(np.array(fsurmargerow).T, columns = list(fsurmargerow.index), index = list(fsurmargerow.columns))
        fmargecol =  self._rdivision(fmargerowT,marge_col)
        Tc =  pd.DataFrame(np.array(fmargecol).T, columns = list(fmargecol.index), index = list(fmargecol.columns)) - 1

        # Weights and svd of Tc
        weightedTc = self._rmultiplication(self._rmultiplication(Tc.T,np.sqrt(marge_col)).T,np.sqrt(marge_row))
        U, s, V = linalg.svd(weightedTc.T, full_matrices=False)
        ncp0 = min(len(weightedTc.iloc[0]),len(weightedTc), self.nf)
        U = U[:,:ncp0]
        V = V.T[:,:ncp0]
        s = s[:ncp0]
        tmp = V
        V = U
        U = tmp
        mult = np.sign(np.sum(V, axis = 0))

        #final V
        mult1 = pd.DataFrame(np.array(pd.DataFrame(np.array(self._rmultiplication(pd.DataFrame(V.T),mult)))).T)
        V = pd.DataFrame()
        for i in range(len(mult1)):
            V[i] = mult1.iloc[i]/np.sqrt(marge_col[i])
        V = np.array(V).T

        #final U
        mult1 = pd.DataFrame(np.array(pd.DataFrame(np.array(self._rmultiplication(pd.DataFrame(U.T),mult)))).T)
        U = pd.DataFrame()
        for i in range(len(mult1)):
            U[i] = mult1.iloc[i]/np.sqrt(marge_row[i])
        U = np.array(U).T
        

        #computing the contribution
        eig = s**2
        for i in range(len(V[0])):
            V[:,i] = V[:,i]*np.sqrt(eig[i])
        coord_col = V

        for i in range(len(U[0])):
            U[:,i] = U[:,i]*np.sqrt(eig[i])
        coord_row = U

        coord_col = coord_col**2

        for i in range(len(coord_col[0])):
            coord_col[:,i] =  (coord_col[:,i]*marge_col)/eig[i]

        contrib =  coord_col*100
        # cols = [str('Dim. ')+str(i) for i in range(ncp0)]
        return contrib