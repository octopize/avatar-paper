from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
import saiph
from numpy.typing import NDArray

from metrics.privacy_metrics.avatars_are_k_hit import avatars_are_k_hit

from metrics.privacy_metrics.hidden_rate import hidden_rate
from metrics.privacy_metrics.local_cloaking import get_local_cloaking


class SecurityMetrics:
    def fit(self, records: pd.DataFrame, avatars: pd.DataFrame, *, nf: int = 5) -> None:
        """Fit the metrics class.

        Arguments:
            records: original data
            avatars: avatarized data
        """
        self._records = records
        self._avatars = avatars
        _, model, param = saiph.fit(records, nf=nf)
        self._model = model
        self._param = param
        self._coord_original = saiph.transform(records, model, param)
        self._coord_avatar = saiph.transform(avatars, model, param)

        are_first_hit = avatars_are_k_hit(
            self._coord_original, self._coord_avatar, distance_metric="minkowski", k=1
        )
        # Hidden rate
        self.hidden_rate = hidden_rate(are_first_hit)

        # Local cloaking
        self.local_cloaking = get_local_cloaking(
            self._coord_original,
            self._coord_avatar,
        )
