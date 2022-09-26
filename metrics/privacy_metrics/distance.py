from typing import Any, List, Optional, Tuple

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from metrics.faiss_knn import FaissKNeighbors

def get_distances_closest_records(original: pd.DataFrame, synthetic: pd.DataFrame, searching_frame: int)->None: 
    nn = FaissKNeighbors(k=searching_frame)
    nn.fit(np.array(synthetic))

    # index.search returns two arrays (distances, indices)
    # https://github.com/facebookresearch/faiss/wiki/Getting-started

    distances, indices = nn.index.search(
        original.to_numpy().astype(np.float32), k=searching_frame
    )

    indices_distances = list(zip(indices, distances))
    return indices_distances