import pandas as pd

from lsg.dimension.projection import Projection
from lsg.security_metrics.avatars_are_first_hit import avatars_are_first_hit
from lsg.security_metrics.hidden_rate import hidden_rate
from lsg.security_metrics.local_cloaking import local_cloaking
from lsg.security_metrics.record_to_avatar_distance import record_to_avatar_distance


class Security_metrics:
    def __init__(self):
        pass

    def fit(self, records_set, avatars_set, nf):
        """Fit the metrics class.

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

    def correlation_protection_rate(self, variable_list):
        """Apply the correlation protection rate metric.

        Arguments:
            variable_list {list} -- name of the variables known by the attacker
        """
        are_first_hit = avatars_are_first_hit(
            self._coord_original, self._coord_avatar, distance_metric="minkowski"
        )
        protection_rate = hidden_rate(are_first_hit)

        return correlation_protection_rate(
            self._records_set,
            self._avatars_set,
            variable_list,
            worst_case_protection_rate=protection_rate,
        )

    def inference_metrics(self, variable_list, target=None):
        """Apply the attack by inference metric.

        Arguments:
            variable_list {list} -- name of the variables known by the attacker
        """
        return inference_metrics(
            self._records_set, self._avatars_set, variable_list, target
        )

    def fit_transform(self, records_set, avatars_set, variable_list):
        """Apply all the metrics and builds a report dataframe.

        Arguments:
            records_set {dataframe} -- original data
            avatars_set {dataframe} -- avatarized data
            variable_list {list} -- name of the variables known by the attacker

        Returns:
            dataframe -- report of the metrics
        """
        # Projection
        pr = Projection()
        coord_original, _ = pr.fit_transform(records_set)
        coord_avatar = pr.transform(avatars_set)

        # Preparation
        are_first_hit = avatars_are_first_hit(
            coord_original, coord_avatar, distance_metric="minkowski"
        )
        distances = record_to_avatar_distance(coord_original, coord_avatar)

        protection_rate = hidden_rate(are_first_hit)

        local_cloak = local_cloaking(
            coord_original, coord_avatar, distances, subset_indices=None
        )

        corr_protection = correlation_protection_rate(
            records_set,
            avatars_set,
            variable_list,
            worst_case_protection_rate=protection_rate,
        )

        # FIXME: target is commented because it is undefined
        inference_attack = inference_metrics(
            records_set,
            avatars_set,
            variable_list,  # target
        )

        d = {
            "rate": [
                protection_rate,
                corr_protection,
                (1 - local_cloak["empty_rate"]) * 100,
                inference_attack["inference_hidden_rate"] * 100,
            ],
            "median": [
                None,
                None,
                (local_cloak["avatars_median"] + local_cloak["records_median"]) / 2,
                inference_attack["inference_local_cloaking"],
            ],
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
