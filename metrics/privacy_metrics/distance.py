import numpy as np
import pandas as pd
from numpy.typing import NDArray

from metrics.faiss_knn import FaissKNeighbors


def get_distances_closest_records(
    records: pd.DataFrame, synthetic: pd.DataFrame, searching_frame: int
) -> List[Tuple[NDArray[np.int_], NDArray[np.float_]]]:
    """Get index and distances of the closest records.

    Arguments
    ---------
        records: Original records
        synthetic: Synthetic data
        searching_frame: number of neighbor to found
    Returns
    -------
        indices_distances: indices, distances nearest neighbor among original records.

    """
    nn = FaissKNeighbors(k=searching_frame)
    nn.fit(np.array(records))

    distances, indices = nn.index.search(
        synthetic.to_numpy().astype(np.float32), k=searching_frame
    )

    indices_distances = list(zip(indices, distances))
    return indices_distances
