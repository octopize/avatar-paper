import pytest
import pandas as pd
import numpy as np
from metrics.privacy_metrics.tmp_dcr_nndr import compare, prepare_common_data_format

from metrics.privacy_metrics.dcr_nndr import get_dcr, get_nndr


original_without_duplicates_1 = pd.DataFrame(
    {
        "variable_1": [1,2,9,10],
        "variable_2": ["Oui", "Oui", "Non", "Non"],
        "id": [1,2,3,4],
    }
)

synthetic_without_duplicates_1 = pd.DataFrame(
    {
        "variable_1": [1,2,9,10],
        "variable_2": ["Oui", "Oui", "Non", "Non"],
        "id": [1,2,3,4],
    }
)

synthetic_without_duplicates_2 = pd.DataFrame(
    {
        "variable_1": [3,4,11,12],
        "variable_2": ["Oui", "Oui", "Non", "Non"],
        "id": [1,2,3,4],
    }
)

original_with_duplicates_2 = pd.DataFrame(
    {
        "variable_1": [1,1,9,9],
        "variable_2": ["Oui", "Oui", "Non", "Non"],
        "id": [1,2,3,4],
    }
)

synthetic_with_duplicates_3 = pd.DataFrame(
    {
        "variable_1": [2,2,10,10],
        "variable_2": ["Oui", "Oui", "Non", "Non"],
        "id": [1,2,3,4],
    }
)

synthetic_with_duplicates_4 = pd.DataFrame(
    {
        "variable_1": [1,1,9,9],
        "variable_2": ["Oui", "Oui", "Non", "Non"],
        "id": [1,2,3,4],
    }
)

synthetic_without_duplicates_5 = pd.DataFrame(
    {
        "variable_1": [3,4,11,12],
        "variable_2": ["Non", "Non", "Oui", "Oui"],
        "id": [1,2,3,4],
    }
)

synthetic_without_duplicates_6 = pd.DataFrame(
    {
        "variable_1": [1,1,9,9],
        "variable_2": ["Non", "Non", "Oui", "Oui"],
        "id": [1,2,3,4],
    }
)

continuous_val = ["variable_1"]
categorical_val = ['variable_2']

@pytest.mark.parametrize(
    "original,synthetic,expected_distance",
    [
        # this first case was verified with a whiteboard
        (original_without_duplicates_1, synthetic_without_duplicates_1, 0),
        (original_without_duplicates_1, synthetic_without_duplicates_2, 0.3222516933177448),
        (original_without_duplicates_1, synthetic_with_duplicates_3, 0),
        (original_without_duplicates_1, synthetic_with_duplicates_4, 0),
        (original_with_duplicates_2, synthetic_without_duplicates_1, 0.10825317547305485),
        (original_with_duplicates_2, synthetic_without_duplicates_2, 0.5412658773652742),
        (original_with_duplicates_2, synthetic_with_duplicates_3, 0.2165063509461097),
        (original_with_duplicates_2, synthetic_with_duplicates_4, 0),
        (original_without_duplicates_1, synthetic_without_duplicates_5, 1.2333775399766171),
        (original_without_duplicates_1, synthetic_without_duplicates_6, 1),
    ],
)
def test_original_synthetic_dcr(original, synthetic, expected_distance) -> None:
    prep_ori = prepare_common_data_format(original,cat_columns= categorical_val,num_columns=continuous_val)
    prep_syn = prepare_common_data_format(synthetic,cat_columns= categorical_val,num_columns=continuous_val)
    privacy = compare(prep_ori, prep_syn, metrics_to_return="privacy-tests")
    assert privacy["DCR"].mean() == expected_distance

@pytest.mark.parametrize(
    "original,synthetic,expected_distance",
    [
        # this first case was verified with a whiteboard
        (original_without_duplicates_1, synthetic_without_duplicates_1, 0),
        (original_without_duplicates_1, synthetic_without_duplicates_2, 0.5833333333333333),
        (original_without_duplicates_1, synthetic_with_duplicates_3, 0),
        (original_without_duplicates_1, synthetic_with_duplicates_4, 0),
        (original_with_duplicates_2, synthetic_without_duplicates_1, 1),
        (original_with_duplicates_2, synthetic_without_duplicates_2, 1),
        (original_with_duplicates_2, synthetic_with_duplicates_3, 1),
        (original_with_duplicates_2, synthetic_with_duplicates_4, 1),
        (original_without_duplicates_1, synthetic_without_duplicates_5, 0.8737211343957376),
        (original_without_duplicates_1, synthetic_without_duplicates_6, 0.823157418648888
),
    ],
)
def test_original_synthetic_nndr(original, synthetic, expected_distance) -> None:
    prep_ori = prepare_common_data_format(original,cat_columns= categorical_val,num_columns=continuous_val)
    prep_syn = prepare_common_data_format(synthetic,cat_columns= categorical_val,num_columns=continuous_val)
    privacy = compare(prep_ori, prep_syn, metrics_to_return="privacy-tests")
    assert privacy["NNDR"].mean() == expected_distance



def test_get_dcr():
    df1 = pd.DataFrame({"variable_1": [1,10]})
    df2 = pd.DataFrame({"variable_1": [1,2]})

    res = get_dcr(df1, df2)
    expected = (np.array([0.]), np.array([64.]))
    assert res == expected


def test_get_nndr():
    df1 = pd.DataFrame({"variable_1": [1,10, 1]})
    df2 = pd.DataFrame({"variable_1": [1, 1, 2]})

    nndr = get_nndr(df1, df2)

    np.testing.assert_almost_equal(nndr, [1, 0.79, 1] , decimal=3)

