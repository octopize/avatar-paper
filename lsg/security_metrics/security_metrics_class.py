import pandas as pd
import saiph

from lsg.security_metrics.avatars_are_first_hit import avatars_are_first_hit
from lsg.security_metrics.hidden_rate import hidden_rate
from lsg.security_metrics.local_cloaking import local_cloaking
from lsg.security_metrics.record_to_avatar_distance import record_to_avatar_distance


class Security_metrics:
    def fit(self, records_set, avatars_set, nf):
        """Fit the metrics class.

        Arguments:
            records_set {dataframe} -- original data
            avatars_set {dataframe} -- avatarized data
        """
        self._records_set = records_set
        self._avatars_set = avatars_set
        _, model, param = saiph.fit(records_set, nf=nf)
        self._coord_original = saiph.transform(records_set, model, param)
        self._coord_avatar = saiph.transform(avatars_set, model, param)

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
        )

    def fit_transform(self, records_set, avatars_set):
        """Apply all the metrics and builds a report dataframe.

        Arguments:
            records_set {dataframe} -- original data
            avatars_set {dataframe} -- avatarized data

        Returns:
            dataframe -- report of the metrics
        """
        # Projection
        _, model, param = saiph.fit(records_set, nf=nf)
        self._coord_avatar = saiph.transform(avatars_set, model, param)
        self._coord_original = saiph.transform(records_set, model, param)

        coord_original = self._coord_original
        coord_avatar = self._coord_avatar 

        # Preparation
        are_first_hit = avatars_are_first_hit(
            coord_original, coord_avatar, distance_metric="minkowski"
        )
        distances = record_to_avatar_distance(coord_original, coord_avatar)

        protection_rate = hidden_rate(are_first_hit)

        local_cloak = local_cloaking(
            coord_original, coord_avatar, distances, subset_indices=None
        )

        d = {
            "rate": [
                protection_rate,
                (1 - local_cloak["empty_rate"]) * 100,
            ],
            "median": [
                None,
                None,
                (local_cloak["avatars_median"] + local_cloak["records_median"]) / 2,
            ]
        }

        data_frame = pd.DataFrame(
            data=d,
            index=[
                "Hidden rate",
                "Correlation protection rate",
                "Local cloaking",
                "Attack by inference",
            ],
        )

        return data_frame
