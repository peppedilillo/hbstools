from pathlib import Path
from typing import NamedTuple

MET = float
# significance, offset
Change = tuple[float, int]
# significance, changepoint_bin, triggertime_bin
Changepoint = tuple[float, int, int]
# significance, changepoint_MET, triggertime_MET
ChangepointMET = tuple[float, MET, MET]


class GTI(NamedTuple):
    """A record for time intervals"""

    start: MET
    end: MET


Dataset = list[tuple[GTI, Path]]


class Event(NamedTuple):
    """A record for transient events"""

    bkg_pre_start: MET
    bkg_pre_end: MET
    start: MET
    end: MET
    bkg_post_start: MET
    bkg_post_end: MET
