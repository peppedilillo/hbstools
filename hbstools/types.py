from typing import NamedTuple

MET = float
Change = tuple[float, int, int]
ChangeMET = tuple[float, MET, MET]
TTI = tuple[float, float, MET, MET, float, float]


class GTI(NamedTuple):
    """A record for holding GTIs"""

    start: MET
    end: MET
