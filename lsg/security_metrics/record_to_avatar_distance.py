import pandas as pd
from numpy import float_
from numpy.typing import NDArray
from sklearn.metrics.pairwise import paired_distances


def record_to_avatar_distance(
    records: pd.DataFrame, avatars: pd.DataFrame
) -> NDArray[float_]:
    """Compute the distance between each record and its avatar.

    Arguments:
        records
        avatars

    Returns:
        distances
    """
    if records.shape[0] != avatars.shape[0]:
        raise ValueError(
            "dimension",
            "Records and avatars dataframes must have the same number of observations",
        )

    # Default is Euclidean distance, we need to implement other distances
    distances: NDArray[float_] = paired_distances(records, avatars, metric="euclidean")

    return distances
