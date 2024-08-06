from pathlib import Path
from typing import NamedTuple

MET = float
# significance, offset
Change = tuple[float, int]
# significance, changepoint_bin, triggertime_bin
Changepoint = tuple[float, int, int]
# significance, changepoint_MET, triggertime_MET
ChangepointMET = tuple[float, MET, MET]
TTI = tuple[MET, MET, MET, MET, MET, MET]


class GTI(NamedTuple):
    """A record for holding time intervals"""

    start: MET
    end: MET


Dataset = list[tuple[GTI, Path]]
