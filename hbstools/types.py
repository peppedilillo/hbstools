from typing import NamedTuple
from pathlib import Path

MET = float
Change = tuple[float, int]
Changepoint = tuple[float, int, int]
ChangepointMET = tuple[float, MET, MET]
TTI = tuple[float, float, MET, MET, float, float]


class GTI(NamedTuple):
    """A record for holding GTIs"""
    start: MET
    end: MET


Dataset = list[tuple[GTI, Path]]
