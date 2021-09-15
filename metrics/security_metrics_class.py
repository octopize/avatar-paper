from avatars_are_first_hit import avatars_are_first_hit
from hidden_rate import hidden_rate
from local_cloaking import local_cloaking
from record_to_avatar_distance import record_to_avatar_distance

from projection import Projection


class Security_metrics:
    def __init__(self):
        pass

    def fit(self, records_set, avatars_set, nf):
        """Fits the metrics class

        Arguments:
            records_set {dataframe} -- original data
            avatars_set {dataframe} -- avatarized data

        """
        self._records_set = records_set
        self._avatars_set = avatars_set
        pr = Projection()
        self._coord_original, _ = pr.fit_transform(records_set, nf=nf)
        self._coord_avatar = pr.transform(avatars_set)

        are_first_hit = avatars_are_first_hit(
            self._coord_original, self._coord_avatar, distance_metric="minkowski"
        )
        # Hidden rate
        self.hidden_rate = hidden_rate(are_first_hit)

        self._distances = record_to_avatar_distance(
            self._coord_original, self._coord_avatar
        )

        # Local cloaking
        self.local_cloaking = local_cloaking(
            self._coord_original,
            self._coord_avatar,
            self._distances,
            subset_indices=None,
        )
