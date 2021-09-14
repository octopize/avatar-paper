from ..dimension.projection import Projection
from ..faiss_knn import FaissKNeighbors
from pandas.api.types import is_numeric_dtype
import numpy as np
from sklearn.metrics import confusion_matrix, accuracy_score, r2_score, mean_squared_error


def inference_metrics(records_set, avatars_set, variable_known, target, k = 1, ncp = None):

    """
    Computes two metrics to evaluate the performance of an attack by inference on an avatarized data set.
    The philosophy of this inference metrics is to evaluate the risk of good prediction using the nearest neighbor of
    each sensitive individual.

    We compute a knn with k=1 to evaluate this risk.
    Outputs metrics are commonly used metrics in machine learning.

    Arguments:
        records_set {dataframe} -- a numpy or pandas dataframe containing numerical, and/or categorical of real values 
        avatars_set {dataframe} -- a numpy or pandas dataframe containing numerical, and/or categorical of avatarized values
        variable_known {list} -- a list containing variables names known by the attacker
        target {string} -- The target variable by the attacker (should be numeric or categorical)
        k {int} -- number of neighbors used for knn, by default k = 1
        ncp {int} -- number of component choosen for data projection, by default ncp = min(5, len(variable_known))
        
    Returns:
        dict {dict} -- A dictonnary containing metrics calculated by the function, depending on target type 
            if target is continious :
                R2_score : coefficient of determination
                RMSE : root-mean-square error
                SMAPE : Symmetric mean absolute percentage error 
            
            if target is categorical :
                acc : accuracy score 
                cf_mat : confusion matrix

        prediction {array} -- The predicted value using knn
        """


    if ncp == None:
        ncp = min(5, len(variable_known))

    if ncp > len(variable_known):
        print('Warning - specifying ncp higher than number of variables known might lead to error if variables are continuous or with few modalities, please consider setting a lower number of ncp')

    if avatars_set.shape[0] != records_set.shape[0]:
        raise ValueError('dimension', 'Records set and avatars set dataframes must have the same number of observations')
    
    if  all(item in records_set.columns.tolist() for item in variable_known) == False : 
        raise ValueError( 'Data and avatar dataframes must have the same columns names')

    if target not in avatars_set.columns:
        raise ValueError('target', 'target must be in records_set.columns')

    if all(item in records_set.columns.tolist() for item in variable_known) == False :
        raise ValueError('variable_known', 'variable_known must be in records_set.columns')

    if target in variable_known : 
        raise ValueError("Your target variable shouldn't be known.")

    if ncp > len(variable_known) : 
        raise ValueError("ncp should be equal or lower to the number of known variables.")
    
    

    # Parameters
    avatar_known = avatars_set[variable_known]
    data_known = records_set[variable_known]
    avatar_target = avatars_set[target]

    # create copy of records_set to work to avoid altering it in return
    working = records_set.copy()

    ### Step 1 : Projection in the known space.

    # fit the space
    pr = Projection()
    avatar_coord, __ = pr.fit_transform(avatar_known, nf = ncp, col_w = None)
    # get indivduals coordinates 
    ori_coord = pr.transform(data_known)
    # get coordinates of nearest avatar(s) neighbor(s) 
    nn = FaissKNeighbors(k = k)
    nn.fit(np.array(avatar_coord))

    # Save index neigbord(s)
    __ , index = nn.predict(np.array(ori_coord))


    ## Step 2 : predict using target type and compute metrics score

    # if target is numeric : 
    if is_numeric_dtype(avatar_target):
        
        # estimate the target value
        prediction = []
        for i in index :
            prediction.append(avatars_set.loc[i, target].mean())
        
        working['estimation'] = prediction
        prediction = np.array(prediction)
        # Compute metrics
        working['symetric_percentage_error'] = abs( working['estimation'] - working[target]) / (abs(working[target]) + abs(working['estimation']) ) * 200

        # to avoid impossible division by 0
        for i in range(working.shape[0]) : 
            # TODO There is a problem here and put conditions in parenthesis doesn't solve it
            if (working.loc[i, 'estimation'] == 0) & (working.loc[i, target] == 0) : 
                working.loc[i, 'symetric_percentage_error']  = 0
                
        # replace infinite value by Nan (when divided by 0)
        working.replace([np.inf, -np.inf], np.nan, inplace=True)
        R2_score = r2_score(working[target], working['estimation'])
        RMSE = mean_squared_error(working[target], working['estimation'])
        SMAPE = working['symetric_percentage_error'].mean(skipna=True)
        metrics = {'R2_score' : R2_score, 'RMSE' : RMSE, 'SMAPE' : SMAPE}


    # if target is categorical : 
    else : 
        # estimate the target value
        avatar_target = avatars_set[target].astype('category').cat.codes
        votes = np.array(avatar_target)[index]
        prediction = np.array([np.argmax(np.bincount(x)) for x in votes])
        working['estimation'] = prediction

        # Compute metrics
        real = working[target].astype('category').cat.codes
        real_numpy = real.values.flatten()

        cf_mat = confusion_matrix(real_numpy, prediction)
        acc = accuracy_score(real_numpy, prediction)

        metrics= {'accuracy' : acc, 'confusion_matrix' : cf_mat}

    return (metrics, prediction)
