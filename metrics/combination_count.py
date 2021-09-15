import numpy as np
import pandas as pd


def combination_count(records_set_data_frame, sensitive_variable):
    """Counts the number of individuals for each combination of sensitive variables

    Arguments:
        records_set_data_frame {dataframe} -- a numpy array given by the function values_counts
        sensitive_variable {array} -- a numpy array giving the index of each sensitive variable in the dataset

    Return:
        dataframe -- a dataframe counting the presence of each combination
    """

    combination = (
        records_set_data_frame.groupby(sensitive_variable)
        .size()
        .reset_index()
        .rename(columns={0: "_count"})
    )

    return combination
