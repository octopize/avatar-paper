import numpy as np
import pandas as pd
from numpy.typing import NDArray, List, Tuple

from metrics.faiss_knn import FaissKNeighbors

def get_distances_closest_records(
    records: pd.DataFrame, synthetic: pd.DataFrame, searching_frame: int
) -> List[Tuple[NDArray[np.int_], NDArray[np.float_]]]:
    """Get index and distances of the closest records.

    Arguments
    ---------
        records: Original records
        synthetic: Synthetic data
        searching_frame: number of neighbors to find
    Returns
    -------
        indices_distances: indices, distances nearest neighbor among original records.

    """
    nn = FaissKNeighbors(k=searching_frame)
    nn.fit(np.array(records))

    # index.search returns two arrays (distances, indices)
    # https://github.com/facebookresearch/faiss/wiki/Getting-started
    distances, indices = nn.predict(synthetic.to_numpy().astype(np.float32))

    indices_distances = list(zip(indices, distances))
    return indices_distances
