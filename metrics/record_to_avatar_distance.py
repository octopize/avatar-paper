from sklearn.metrics.pairwise import paired_distances
import pandas as pd


# def record_to_avatar_distance(full_bind_coordinates, nf = None):
def record_to_avatar_distance(records_set_coordinates, avatars_set_coordinates):
    """ Compute the distance between each record and its avatar
    
    Arguments:
        records_set_coordinates {dataframe} -- a numpy or pandas dataframe containing coordinates of records
        avatars_set_coordinates {dataframe} -- a numpy or pandas dataframe containing coordinates of avatars
    
    Returns:
        dataframe -- a pandas dataframe gathering the distance between each record and its avatar
    """

    if records_set_coordinates.shape[0] != avatars_set_coordinates.shape[0]:
        raise ValueError('dimension', 'Records set and avatars set dataframes must have the same number of observations')

    ## Compute distance
    # Default is Euclidean distance, we need to implement other distances
    distances = paired_distances(records_set_coordinates, avatars_set_coordinates, metric='euclidean')

    return distances
