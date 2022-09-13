from numpy import bool_
from numpy.typing import NDArray


def hidden_rate(are_first_hit: NDArray[bool_]) -> float:
    """Test for each record if the nearest avatar is the one generated by the original record itself.

    Arguments:
        are_first_hit: whether the nearest neighbor of the supplied
          original record is the avatar created from this record

    Returns:
        float: the percentage of records that could be considered "safe"
    """
    res: float = (1 - (are_first_hit.sum() / len(are_first_hit))) * 100
    return res