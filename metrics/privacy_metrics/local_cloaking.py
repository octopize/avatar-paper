from typing import Any, List, NamedTuple, Tuple

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from metrics.faiss_knn import FaissKNeighbors


class LocalCloakingResults(NamedTuple):
    """Store the results of local cloaking metric for intuitive access.

    Returns:
        empty_rate: the percentage of individuals who have no record nor avatar in radius
        avatars_median: median number of avatars between a record and its avatar
        records_median: median number of records between a record and its avatar
        hit_counts: avatar and records count in radius for each individual
    """

    empty_rate: float
    avatars_median: float
    records_median: float
    hit_counts: List[Tuple[int, int]]


def get_local_cloaking(
    records_coordinates: pd.DataFrame,
    avatars_coordinates: pd.DataFrame,
) -> LocalCloakingResults:
    """Get summary of global radius knn metric."""
    if records_coordinates.shape[0] != avatars_coordinates.shape[0]:
        raise ValueError(
            "dimension",
            "Records set and avatars set dataframes must have the same number of observations",
        )

    # High-level approach:
    # For each generated avatar, we get its l closest neighbor
    # For each generated avatar, we check how many datapoints separates it from
    # the original data

    n = records_coordinates.shape[0]

    # We use a single df with original data followed by avatars
    # Avatars have indices ≥ n
    both: NDArray[Any] = np.ascontiguousarray(
        pd.concat([records_coordinates, avatars_coordinates])
        .reset_index(drop=True)
        .to_numpy()
    )

    nn = FaissKNeighbors(k=4000)
    nn.fit(np.array(both))

    # index.search returns two arrays (distances, indices)
    # https://github.com/facebookresearch/faiss/wiki/Getting-started
    distances, indices = nn.index.search(
        records_coordinates.to_numpy().astype(np.float32), k=4000
    )

    # indices will look like this for record index k and n = 10
    # [k, 12, 13, 2, k + n, ...]  # noqa: E800
    # _^----------------------- the closest to k is usually k
    # ____^-------------------- 12 is an avatar because 12 >= n
    # ____________^------------ 2 is a record because 2 < n
    # _______________^--------- k+n-1 is the avatar generated by k

    indices_distances = list(zip(indices, distances))

    # hit_counts is an array of arrays [n_closest_avatars, n_closest_records]
    # separating each original record from its closest avatar
    hit_counts = [
        get_counts(current, indice_distance[0], indice_distance[1], n)
        for current, indice_distance in enumerate(indices_distances)
    ]

    no_hit: NDArray[np.int_] = np.array([sum(t) == 0 for t in hit_counts])

    results = LocalCloakingResults(
        np.mean(no_hit),
        float(np.median([t[0] for t in hit_counts])),
        float(np.median([t[1] for t in hit_counts])),
        hit_counts,
    )

    return results


def get_counts(
    current: int,
    indices: NDArray[np.int_],
    distances: NDArray[np.float_],
    number_individual: int,
) -> Tuple[int, int]:
    avatars_count = 0
    records_count = 0
    duplicated = False

    for k, i in enumerate(indices):
        if i == current:
            # We don't count the current indice (logically)
            # Note that it is not guaranteed to be the first, because in
            # some cases a generated avatar is exactly equal to a
            # previously existing datapoint.
            continue

        # This is the avatar,
        #    -> if the next neighbor distance is equal to the avatar's distance
        #    then continue and look to distance
        if i == current + number_individual or i == -1:
            if i == -1 or k + 1 == len(indices):
                break
            if distances[k] == distances[k + 1]:
                duplicated = True
            else:
                break

        # We don't count his avatar (logically)
        if i == current + number_individual:
            continue
        elif i >= number_individual:
            avatars_count += 1
        else:
            records_count += 1

        # if we pass the avatar's index, we compare the distances (after counting)
        if duplicated:
            if i == -1 or k + 1 == len(indices):
                break
            if distances[k] < distances[k + 1]:
                break

    r = (avatars_count, records_count)

    return r