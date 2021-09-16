def combination_count(records_set_data_frame, sensitive_variable):
    """Count the number of individuals for each combination of sensitive variables.

    Arguments:
        records_set_data_frame {dataframe}
        sensitive_variable {array}

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
