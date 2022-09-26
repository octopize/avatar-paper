from typing import List

import pandas as pd

from metrics.privacy_metrics.distance import get_distances_closest_records


def get_dcr(
    original_coordinates: pd.DataFrame, avatar_coordinates: pd.DataFrame
) -> List[float]:
    """Get distances to the closest records."""
    indices_distances = get_distances_closest_records(
        original_coordinates, avatar_coordinates, searching_frame=1
    )
    _, distances = zip(*indices_distances)
    return [distance[0] for distance in distances]


def get_nndr(
    original_coordinates: pd.DataFrame, avatar_coordinates: pd.DataFrame
) -> List[float]:
    """Get nearest neighbors distance ratio."""
    indices_distances = get_distances_closest_records(
        original_coordinates, avatar_coordinates, searching_frame=2
    )
    _, distances = zip(*indices_distances)

    ratio = [
        1 if distance[1] == 0 else distance[0] / distance[1] for distance in distances
    ]
    return ratio
