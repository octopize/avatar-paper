from  metrics.privacy_metrics.distance import get_distances_closest_records


def get_dcr(original_coordinates, avatar_coordinates):
    indices_distances = get_distances_closest_records(original_coordinates, avatar_coordinates, searching_frame=1)
    _, distances = zip(*indices_distances)
    return distances

def get_nndr(original_coordinates, avatar_coordinates):
    indices_distances = get_distances_closest_records(original_coordinates, avatar_coordinates, searching_frame=2)
    _, distances = zip(*indices_distances)

    ratio = [1 if distance[1] == 0 else distance[0]/distance[1] for distance in distances]
    return ratio
