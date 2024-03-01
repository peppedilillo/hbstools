from functools import wraps

import pandas as pd

from hbstools.types import ChangepointMET
from hbstools.types import GTI
from hbstools.types import MET
from hbstools.types import TTI

FCOLS = [
    "bkg_pre_start",
    "bkg_pre_end",
    "event_start",
    "event_end",
    "bkg_post_start",
    "bkg_post_end",
]


def compute_bkg_pre(
    trigtime: MET,
    changepoint: MET,
    gti: GTI,
    ends_seconds: float,
    duration_seconds: float,
) -> tuple[float, float]:
    """Times from changepoint of an interval ending `m` steps behind the triggertime,
    with duration binning / alpha."""
    if trigtime - ends_seconds - duration_seconds < gti.start:
        start = -(changepoint - gti.start)
        end = -ends_seconds + (trigtime - changepoint)
    else:
        end = (trigtime - changepoint) - ends_seconds
        start = end - duration_seconds
    return start, end


def compute_bkg_post(
    trigtime: MET,
    changepoint: MET,
    gti: GTI,
    start_seconds: float,
    duration_seconds: float,
) -> tuple[float, float]:
    """Times from changepoint of an interval starting after skip steps from triggertime,
    with duration binning / alpha."""
    if trigtime + start_seconds + duration_seconds > gti.end:
        end = gti.end - changepoint
        start = max(end - duration_seconds, trigtime - changepoint)
    else:
        start = (trigtime - changepoint) + start_seconds
        end = start + duration_seconds
    return start, end


def _format(
    preinterval_duration_seconds: float,  # (binning / alpha)
    preinterva_ends_seconds: float,  # (binning * m)
    postinterval_duration_seconds: float,  # binning / alpha
    postinterval_start_seconds: float,  # binning * skip
):
    def _format_result(
        result: ChangepointMET,
        gti: GTI,
    ) -> TTI:
        """Transforms a single focus results (times expressed as mets) into a TTI event formatted like:
        (bkg_pres_start, bkg_pre_ends, event_starts, event_ends, bkg_post_start, bkg_post_end)
        """
        significance, changepoint, trigtime = result
        bkg_pre = compute_bkg_pre(
            trigtime,
            changepoint,
            gti,
            preinterva_ends_seconds,
            preinterval_duration_seconds,
        )
        bkg_post_start, bkg_post_end = compute_bkg_post(
            trigtime,
            changepoint,
            gti,
            postinterval_start_seconds,
            postinterval_duration_seconds,
        )
        event_interval = changepoint, changepoint + bkg_post_start
        # noinspection PyTypeChecker
        return *bkg_pre, *event_interval, bkg_post_start, bkg_post_end

    return _format_result


def as_dataframe(func):
    @wraps(func)
    def wrapper(*args, **kwargs) -> pd.DataFrame:
        return pd.DataFrame(
            func(*args, **kwargs),
            columns=FCOLS,
        )

    return wrapper


@as_dataframe
def format_results(
    results: dict[GTI, list[ChangepointMET]],
    intervals_duration_seconds: float,  # (binning / alpha)
    preinterva_ends_seconds: float,  # (binning * m)
    postinterval_start_seconds: float,  # binning * skip
) -> list[TTI]:
    format_result = _format(
        intervals_duration_seconds,
        preinterva_ends_seconds,
        intervals_duration_seconds,
        postinterval_start_seconds,
    )
    return [
        format_result(changepoint, gti)
        for gti, cps in results.items()
        for changepoint in cps
    ]
