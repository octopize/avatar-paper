import pandas as pd
import numpy as np

from .combination_count import combination_count


def correlation_protection_rate(
    records_set,
    avatars_set,
    variable_list=None,
    worst_case_protection_rate=None,
    boundaries=0.05,
):
    """Calculate the protection rate regarding values held by variables likely to be found by an attacker

    Arguments:
        records_set {dataframe} -- a pandas dataframe containing numerical, and/or categorical values
        avatars_set {dataframe} -- a pandas dataframe containing numerical, and/or categorical values
        variable_list {array} -- list containing the index of the variables you want to check
        worst_case_protection_rate {float} -- a float returned by the hidden_rate function
        boundaries {float} -- a float giving percentage of variation for which we consider two values as identical

    Raises:
        ValueError: float within [0,100] needs to be specified

    Returns:
        float -- the proportion of records that could be considered "safe" in a probable scenario
    """

    if worst_case_protection_rate is None:
        raise ValueError(
            "worst_case_protection_rate", "float within [0,100] needs to be specified"
        )

    # Preprocessing
    working_df = records_set[variable_list]
    N = records_set.shape[0]

    # Count every combination of selected variables in original dataset
    combination_records = combination_count(
        records_set, sensitive_variable=variable_list
    )

    # keep only unique combinations
    unique_records = combination_records[combination_records["_count"] == 1]

    n = unique_records.shape[0]
    v = working_df.shape[1]
    number_matching = []

    # for each unique original record, find those with only one approaching avatar record

    # n = number of unique
    # n2 = number of unique that have a unique similar avatar
    # n3 = n2 * (1-protection rate)

    for i in range(0, n):
        iteration_dataframe = avatars_set[variable_list]
        for c in range(0, v):
            if (
                unique_records.iloc[:, c].dtype == np.float64
                or unique_records.iloc[:, c].dtype == np.int64
            ):
                iteration_dataframe = iteration_dataframe.loc[
                    (
                        iteration_dataframe.iloc[:, c]
                        > (unique_records.iloc[i, c])
                        - (working_df.iloc[:, c].mean() * (boundaries))
                    )
                    & (
                        iteration_dataframe.iloc[:, c]
                        < (unique_records.iloc[i, c])
                        + (working_df.iloc[:, c].mean() * (boundaries))
                    )
                ]
            else:
                iteration_dataframe = iteration_dataframe.loc[
                    iteration_dataframe.iloc[:, c] == unique_records.iloc[i, c]
                ]
        number_matching.append(iteration_dataframe.shape[0])

    n2 = number_matching.count(1)

    # Percentage of protection
    percentage = (1 - round(n2 * (1 - (worst_case_protection_rate / 100))) / N) * 100
    return percentage
