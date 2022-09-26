from typing import Any, List, Optional, Tuple

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from metrics.faiss_knn import FaissKNeighbors


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


