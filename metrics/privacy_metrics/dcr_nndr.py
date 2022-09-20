"""
Provide methods that calculate the distance between two (encoded) sequential datasets:
Privacy Metrics:
* distance to closest records - quantiles
* nearest neighbor distance ratio - quantiles
"""

import pandas as pd
from pandas.api.types import is_numeric_dtype,is_categorical_dtype
import numpy as np
import random
from itertools import product
from pandas import DataFrame, Series
from typing import List, Tuple, Dict, Callable, Any
import scipy.stats as ss
from sklearn.neighbors import NearestNeighbors

def get_categorical_continuous(data: pd.DataFrame):
    categorical_val = []
    continous_val = []
    for column in data.columns:
        if len(data[column].unique()) <= 10:
            categorical_val.append(column)
        else:
            continous_val.append(column)
    
    return categorical_val, continous_val

def prepare_common_data_format(
    fp: Any, cat_columns: List[str] = [], num_columns: List[str] = []
) -> DataFrame:
    """
    Cast cat columns into pd.Categorical. Cast num columns into numeric. Can accept parquet + csv files
    Any columns not passed in either list will be automatically casted. Anything not numeric is category.
    Set index to be id and sequence pos
    :param fp: filepath or pandas dataframe
    :param cat_columns: List of category columns
    :param num_columns: List of numeric columns
    :returns: same Dataframe with types
    """
    df = fp

    columns = list(df.columns)
    input_columns = cat_columns + num_columns

    assert "id" in columns, "id col not in dataframe. please rename if possible"

    columns.sort()
    input_columns.sort()

    if input_columns == columns:
        cat_columns_convert = cat_columns
        num_columns_convert = num_columns
    else:
        unlabeled_columns = [col for col in columns if col not in input_columns]
        cat_columns_convert = cat_columns
        num_columns_convert = num_columns

        for unlabeled_column in unlabeled_columns:
            if is_numeric_dtype(df[unlabeled_column]):
                num_columns_convert.append(unlabeled_column)
            else:
                cat_columns_convert.append(unlabeled_column)

    for col in columns:
        if col != "id":
            if col in cat_columns_convert:
                df.loc[:, col] = pd.Categorical(df[col])
                # nan needs to be considered a category
                if df.loc[:, col].isnull().sum() != 0:
                    df.loc[:, col] = df.loc[:, col].cat.add_categories("").fillna("")
            elif col in num_columns_convert:
                # throw error when nan in cols
                assert (
                    df.loc[:, col].isnull().sum() == 0
                ), "Numeric columns contain NaN values"
                df.loc[:, col] = pd.to_numeric(df[col])

    df["sequence_pos"] = df.groupby("id").cumcount()

    return df.set_index(["id", "sequence_pos"])


def mixed_distance(x,y,categoric_slice):
    n = x.shape[0]
    i = 0
    res = 0
    for i in range(categoric_slice):
        res += abs(x[i] != y[i])
    for i in range(categoric_slice, n):
        res += abs(x[i] - y[i])
    return res
    
def _generate_column_type_dictionary(data) -> dict:
    column_types = dict(zip(data.dtypes.index, data.dtypes))
    column_type_dict = {}
    true_vector = []
    for col_name, dtype in column_types.items():
        if is_numeric_dtype(dtype):
            column_type_dict[col_name] = "number"
            true_vector.append(True)
        elif is_categorical_dtype(dtype):
            column_type_dict[col_name] = "category"
            true_vector.append(True)
        else:
            true_vector.append(False)
    assert sum(true_vector) == len(
        data.columns
    ), " Data must only have number and category column types"
    return column_type_dict

def check_common_data_format(data):
    """
    Perform validation for common data format. Used to screen input to validate if data is in common data format
    Requirements:
    * Pandas DataFrame
    * numerical and cat dtypes
    * must have id and column as multi-index
    :param data: Pandas Dataframe
    :raise: AssertionError if any criteria is not met
    """

    assert isinstance(data, DataFrame), "Data is not Pandas DataFrame"

    column_types = dict(zip(data.dtypes.index, data.dtypes))

    column_type_dict = {}
    true_vector = []

    for col_name, dtype in column_types.items():
        if is_numeric_dtype(dtype):
            column_type_dict[col_name] = 'numeric'
            true_vector.append(True)
        elif is_categorical_dtype(dtype):
            column_type_dict[col_name] = 'categorical'
            true_vector.append(True)
        else:
            true_vector.append(False)

    assert sum(true_vector) == len(data.columns), "Data must only have numeric and categorical column types"

    index_names = data.index.names
    assert "id" in index_names, "id not in Index"


def _bin_data(target_col, syn_col, col_type: str, number_of_bins: int):
    """
    Bin single target/synthetic column.
    """
    if col_type == "number":
        cardinality = target_col.nunique()
        if cardinality > number_of_bins:
            filled_bins = 0
            q = number_of_bins
            while filled_bins < number_of_bins:
                binned_target, bins = pd.qcut(
                    target_col, q=q, retbins=True, duplicates="drop", precision=1
                )
                filled_bins = sum(binned_target.value_counts() != 0)
                q += 1
        else:
            bins = sorted(target_col.unique())
        binned_target = pd.cut(target_col, bins=bins, include_lowest=True, precision=1)
        # if synthetic has categories not in target, na value will appear, that we drop
        binned_syn = pd.cut(
            syn_col, bins=bins, include_lowest=True, precision=1
        ).dropna()
    elif col_type == "category":
        cardinality = target_col.nunique()
        if cardinality > number_of_bins:
            top_cat = (
                target_col.value_counts()
                .sort_values(ascending=False)[: number_of_bins - 1]
                .index.to_list()
            )
            other = (
                target_col.value_counts()
                .sort_values(ascending=False)[number_of_bins - 1 :]
                .index.to_list()
            )
            first_dictionary = {cat: cat for cat in target_col.unique().to_list()}
            first_dictionary.update({cat: "*" for cat in other})
            binned_target = target_col.map(first_dictionary)
            # if synthetic has categories not in target, na value will appear, that we drop
            binned_syn = syn_col.map(first_dictionary).dropna()
        else:
            binned_target = target_col
            binned_syn = syn_col
    else:
        raise Exception("Col type not recognized")
    return binned_target, binned_syn


def _bin_looped(target, synthetic, column_dictionary, number_of_bins):
    t = pd.DataFrame()
    s = pd.DataFrame()
    for col, col_type in column_dictionary.items():
        binned_target, binned_syn = _bin_data(
            target[col], synthetic[col], col_type, number_of_bins
        )
        t[col] = binned_target
        s[col] = binned_syn
    return t, s


def _flatten_table(data: DataFrame, column_type_dictionary: dict) -> DataFrame:
    """
    Flatten a table from long format into wide format.
    Long format
    | id | col_1 | col_2 | record_pos | sequence_pos |
    |----|-------|-------|------------|--------------|
    | 1  | a     | 10    | 1          | 3            |
    | 1  | b     | 20    | 2          | 4            |
    | 2  | c     | 30    | 1          | 0            |
    | 2  | d     | 40    | 2          | 1            |
    Wide format
    | id | col_a_1 | col_a_2 | col_b_1 | col_b_2 |
    |----|---------|---------|---------|---------|
    | 1  | a       | b       | 10      | 20      |
    | 2  | c       | d       | 30      | 40      |
    This is used as a utility function to do data preprocessing for metric calculations.
    Note: Sequence pos and record pos is dropped from returned dataframe.
    :param data: Pandas Dataframe
    :param column_type_dictionary: dict mapping columns to types. Note, pandas.pivot() coerces types so we need to reassign them.
    :returns: wide Pandas DataFrame
    """

    to_exclude = ["sequence_pos", "record_pos", "id"]
    col_to_melt = [x for x in data.columns if x not in to_exclude]
    pivot = data.pivot(index="id", columns="sequence_pos", values=col_to_melt)

    # Lots of extra steps because pandas.pivot does not preserve column type

    new_columns = ["{}_{}".format(x[0], x[1]) for x in pivot.columns]
    original_columns_family = ["{}".format(x[0]) for x in pivot.columns]
    new_columns_to_family = dict(zip(new_columns, original_columns_family))
    pivot.columns = new_columns

    for col in pivot:
        original_col_type = column_type_dictionary[new_columns_to_family[col]]
        if original_col_type == "category":
            pivot[col] = pd.Categorical(pivot[col])
        elif original_col_type == "number":
            pivot[col] = pd.to_numeric(pivot[col])

    return pivot


def _prepare_data_for_privacy_metrics(
    tgt_data: DataFrame,
    syn_data: DataFrame,
    column_dictionary: Dict,
    smoothing_factor: float,
) -> Tuple[DataFrame, DataFrame]:
    """
    Data preparation for privacy metrics
    For categorical, ordinal encoding based on joint set of target data and synthetic data.
    For numeric encoding, missing value are imputed with mean and standardized
    :param tgt_data: pandas dataframe
    :param syn_data: pandas dataframe
    :param column_dictionary: column to type mapping
    :param smoothing_factor: smoothing factor
    :returns: privacy ready target + synthetic  (pamdas DataFrame)
    """

    tgt_data_p = tgt_data.copy(deep=True)
    syn_data_p = syn_data.copy(deep=True)

    for column_name, column_type in column_dictionary.items():
        if column_type == "category":

            tgt_data_p[column_name] = tgt_data_p[column_name].cat.codes
            syn_data_p[column_name] = syn_data_p[column_name].cat.codes

        elif column_type == "number":
            # fill na data with mean
            tgt_data_p[column_name] = tgt_data_p[column_name].fillna(
                tgt_data_p[column_name].dropna().mean()
            )
            syn_data_p[column_name] = syn_data_p[column_name].fillna(
                syn_data_p[column_name].dropna().mean()
            )

            # standardize
            syn_data_p[column_name] = (
                syn_data_p[column_name] - tgt_data_p[column_name].mean()
            ) / np.max([tgt_data_p[column_name].std(), smoothing_factor])
            tgt_data_p[column_name] = (
                tgt_data_p[column_name] - tgt_data_p[column_name].mean()
            ) / np.max([tgt_data_p[column_name].std(), smoothing_factor])
        else:
            raise Exception(f"{column_type} Type not supported")

    # drop id col since it's not needed
    tgt_data_p = tgt_data_p.reset_index().drop("id", 1)
    syn_data_p = syn_data_p.reset_index().drop("id", 1)

    return tgt_data_p, syn_data_p


def _get_nn_model(train: DataFrame, cat_slice: int) -> Tuple[np.ndarray]:
    """
    Find nearest neighbors of test in train with first categoric_slice-many variables being categorical.
    :param train: train pandas dataframe
    :param cat_slice: where do category columns end
    :returns: scikit learn nearest_neighbor_model fit with train data
    """
    nearest_neighbor_model = NearestNeighbors(
        metric=lambda x, y: mixed_distance(x, y, cat_slice),
        algorithm="ball_tree",
        n_jobs=None,
    )
    nearest_neighbor_model.fit(train)

    return nearest_neighbor_model


def _calculate_dcr_nndr(
    tgt_data: DataFrame,
    syn_data: DataFrame,
    column_dictionary: dict,
    smoothing_factor: int,
) -> Tuple[Dict, Dict]:
    """
    Function to calculate dcr and nndr. Since DCR and NNDR are related, both are calculated at the same time.
    :param tgt_data: privacy prepared dataset target
    :param syn_data: privacy prepared dataset synthetic
    :param column_dictionary: column to type mapping
    :param smoothing_factor: smoothing factor to avoid small division
    :returns: bins and histograms of dcr and nndr
    """

    # derive model based on sample of features
    model_columns = list(column_dictionary.keys())
    feature_columns = model_columns

    # shift columns to put category features first for distance metric
    category_columns = [
        x for x in feature_columns if column_dictionary[x] == "category"
    ]
    ordered_columns = category_columns + [
        x for x in feature_columns if column_dictionary[x] == "number"
    ]
    # where do cat columns begin
    cat_slice = len(category_columns)

    tgt_data = tgt_data[ordered_columns]
    syn_data = syn_data[ordered_columns]

    assert all(
        tgt_data.columns == tgt_data.columns
    ), "Train and Syn have mismatched columns"

    # training model
    nn_model = _get_nn_model(tgt_data, cat_slice)

    tgt_query_neighbors = nn_model.kneighbors(tgt_data, n_neighbors=2)
    syn_query_neightbors = nn_model.kneighbors(syn_data, n_neighbors=2)

    # Calculating DCR NNDR
    query_dict = {"tgt": tgt_query_neighbors, "syn": syn_query_neightbors}

    privacy_data = {}

    for label, query in query_dict.items():
        if label == 'syn':
            dcr = query[0][:, 0]
            nndr = query[0][:, 0] / np.maximum(query[0][:, 1], smoothing_factor)
            mask = query[0][:, 1]!=0
            nndr[~mask]=1
            df_privacy = pd.DataFrame({"DCR": dcr, "NNDR": nndr})
            privacy_data[label] = df_privacy

    return df_privacy

    #return checks, privacy_tests, bootstrap_df

def _calculate_privacy_tests(tgt_data: DataFrame, syn_data: DataFrame):
    """
    Compute privacy tests for a given target and synthetic data set
    """

    assert sorted(tgt_data.columns) == sorted(
        syn_data.columns
    ), "Target and Synthetic have different columns"

    tgt_dict = _generate_column_type_dictionary(tgt_data)
    syn_dict = _generate_column_type_dictionary(syn_data)
    assert tgt_dict == syn_dict, "Target and Synthetic have different types"

    flat_table_target = _flatten_table(tgt_data.reset_index(), tgt_dict)
    flat_table_syn = _flatten_table(syn_data.reset_index(), tgt_dict)

    # columns now include 1st and 2nd record
    column_dict = _generate_column_type_dictionary(flat_table_target)

    smoothing_factor = 1e-8
    tgt_data_p, syn_data_p = _prepare_data_for_privacy_metrics(
        flat_table_target, flat_table_syn, column_dict, smoothing_factor
    )

    df_privacy = _calculate_dcr_nndr(
        tgt_data_p, syn_data_p, column_dict, smoothing_factor
    )
    return df_privacy


def compare(
    tgt_data: DataFrame,
    syn_data: DataFrame,
    metrics_to_return: List = ["statistical-distances", "privacy-tests"],
) -> Dict:
    """
    Compare a target and synthetic dataset and returns metrics
    :param tgt_data: target data
    :param syn_data: synthetic data
    :param metrics_to_return: list of metric types to return, if empty compute all
    ['statistical-distances','privacy-tests']
    """
    check_common_data_format(tgt_data)
    check_common_data_format(syn_data)

    if "privacy-tests" in metrics_to_return:
        privacy_df = _calculate_privacy_tests(tgt_data, syn_data)

    return privacy_df