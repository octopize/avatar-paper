from typing import Any, List, Optional, Tuple

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from lsg.faiss_knn import FaissKNeighbors


def get_nearest_neighbor(
    original_coordinates: pd.DataFrame,
    avatars_coordinates: pd.DataFrame,
    searching_frame: Optional[int] = None,
) -> List[Tuple[NDArray[np.int_], NDArray[np.float_]]]:
    """Compute nearest neighbor distances and index in both original and avatars data.

    Arguments
    ---------
        original: original coordinates after dimension reduction.
        avatars: avatars coordinates after dimension reduction.
        searching_frame: size of the frame for KNN research.

    Returns
    -------
        indices_distances: indices, distances nearest neighbor among original and avatars.
    """
    if original_coordinates.shape[0] != avatars_coordinates.shape[0]:
        raise ValueError(
            "dimension",
            "original set and avatars set dataframes must have the same number of observations",
        )

    if not searching_frame:
        searching_frame = min(100, len(original_coordinates) * 2)

    # For each generated avatar, we get its l closest neighbor
    # among original and avatars data

    original_coordinates.shape[0]

    # We use a single df with original data followed by avatars
    # Avatars have indices ≥ n
    both: NDArray[Any] = np.ascontiguousarray(
        pd.concat([original_coordinates, avatars_coordinates])
        .reset_index(drop=True)
        .to_numpy()
    )

    nn = FaissKNeighbors(k=searching_frame)
    nn.fit(np.array(both))

    # index.search returns two arrays (distances, indices)
    # https://github.com/facebookresearch/faiss/wiki/Getting-started

    distances, indices = nn.index.search(
        original_coordinates.to_numpy().astype(np.float32), k=searching_frame
    )

    indices_distances = list(zip(indices, distances))

    return indices_distances


def get_closeness_value(
    record: int, indice_distance: Tuple[NDArray[np.int_], NDArray[np.float_]], n: int
) -> int:
    """Get closeness of generated avatar for a specific record.

    Arguments
    ---------
        record: index of the record.
        indice_distance: index and distance of nearest neighbor of a record
        n: overall number of records

    Returns
    -------
        closeness: position of generated avatar among other.
    """
    # Keep only avatars for hit count
    indices = [i for i in indice_distance[0] if i >= n]
    # Get index of avatars value
    avatars_index = [
        x for x in range(len(indice_distance[1])) if indice_distance[0][x] >= n
    ]
    # Use avatars_index to keep only avatars distances
    distances = [indice_distance[1][x] for x in avatars_index]
    # Get position index of generated avatar
    position = [x for x in range(len(indices)) if indices[x] == record + n]
    # If avatar not in list set max value (size of searching_frame)
    if not position:
        value = len(indices)
    else:
        value = position[0]
        # If the generated avatar is at an equal distance to other avatars
        # These are included in radius
        while value < (len(distances) - 1) and distances[value] == distances[value + 1]:
            value += 1
    return value


def get_avatar_closeness(
    indices_distances: List[Tuple[NDArray[np.int_], NDArray[np.float_]]],
) -> List[int]:
    """Get closeness of generated avatar for each record.

    Arguments
    ---------
        indices_distances: index and distances of nearest neighbor for each record.

    Returns
    -------
        closeness: list of generated avatar positions among other for each record.
    """
    n = len(indices_distances)
    closeness = [
        get_closeness_value(current, indice_distance, n=n)
        for current, indice_distance in enumerate(indices_distances)
    ]
    return closeness


def get_unique_index(
    df: pd.DataFrame,
) -> List[int]:
    """Get index of unique records only.

    Arguments
    ---------
        df: dataframe used for unique selection

    Returns
    -------
        list: position unique records in dataframe index.
    """
    df = df.reset_index(drop=True)
    df = df.drop_duplicates(keep=False)
    return list(df.index)
