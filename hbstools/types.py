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
    stop: MET


# Dataset are nothing special really, they are lists of tuples built like this:
# ┌──────────────────────────────────────────────┐┌──────────────────────────────┐
# │[((Path('file1.evt'), Path('file1_gti.fits')),││ GTI(start=0.0, end=54.0)),   │
# │ ((Path('file2.evt'), Path('file2_gti.fits')),││ GTI(start=51.0, end=79.0)),  │
# │ ((Path('file2.evt'), Path('file2_gti.fits')),││ GTI(start=83.0, end=108.0)), │
# │ ((Path('file3.evt'), Path('file3_gti.fits')),││ GTI(start=108.5, end=133.0))]│
# └──────────────────datafiles───────────────────┘└──────────────gtis────────────┘
Dataset = list[tuple[tuple, GTI]]


class Event(NamedTuple):
    """A record for transient events"""

    bkg_pre_start: MET
    bkg_pre_stop: MET
    start: MET
    stop: MET
    bkg_post_start: MET
    bkg_post_stop: MET
