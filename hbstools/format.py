"""
# EVENT FORMATTING RULES

This is a schematic representing the `bkg_pre`, `event` and `bkg_post` intervals.


  BKG_PRE_INTERVAL                                             BKG_POST_INTERVAL
(*<---pre_delta--->*)|....pre_t....|.........post_t.........|(*<---post_delta--->*)|
                              trigger_time                  |
                                   |                        |
                      changepoint  |                        |
                           |(*<-----post_t + tt - cp----->*)|
                                    EVENT_INTERVAL



Given these durations, the interval start and end times are as follow:
    * bkg_pre_start = trigger_time - (pre_t + pre_delta)
    * bkg_pre_end = trigger_time - pre_t
    * event_start = changepoint
    * event_end = trigger_time + post_t
    * bkg_post_start = trigger_time + post_t
    * bkg_post_end = trigger_time + (post_t + post_delta)

For a single exponential smoothing algorithm the interval timings are chosen so that:
    * pre_delta = post_delta = (binning / alpha)
    * pre_t = binning * m
    * post_t = binning * skip

Applying to the definition just above we have the following interval time, as a function
of the algorithm's parameters:
    * bkg_pre_start = trigger_time - binning * (m + 1 / alpha)
    * bkg_pre_end = trigger_time - binning * m
    * event_start = changepoint
    * event_end = trigger_time + binning * skip
    * bkg_post_start = trigger_time + binning * skip
    * bkg_post_end = trigger_time + binning * (skip + 1/alpha)

## Caveats:
    1. it is possible for the changepoint to predate the end of the bkg interval if
       an algorithms is launched with paramters `t_max` larger than `m`.
    2. it is guaranteed that all intervals will fit in the event's GTI. if an interval
       overcome the GTI boundaries, it is "squashed" enough to fit.
    3. if an event occurs across different datafiles, and the GTI of the datafiles are
       adjacent (overlap), the interval is not squashed.
"""

from typing import Callable

from hbstools.types import GTI, MET, Event, ChangepointMET


def compute_bkg_pre(
    trigtime: MET,
    gti: GTI,
    t: float,
    delta: float,
) -> tuple[MET, MET]:
    """Times from changepoint of an interval ending `m` steps behind the triggertime,
    with duration binning / alpha."""
    if trigtime - t - delta < gti.start:
        start = gti.start
        end = - t + trigtime
    else:
        end = trigtime - t
        start = end - delta
    return start, end


def compute_bkg_post(
    trigtime: MET,
    gti: GTI,
    t: float,
    delta: float,
) -> tuple[MET, MET]:
    """Times from changepoint of an interval starting after skip steps from triggertime,
    with duration binning / alpha."""
    if trigtime + t + delta > gti.end:
        end = gti.end
        start = max(end - delta, trigtime)
    else:
        start = trigtime + t
        end = start + delta
    return start, end


def _format(
    pre_delta: float,  # (binning / alpha)
    pre_t: float,  # (binning * m)
    post_delta: float,  # binning / alpha
    post_t: float,  # binning * skip
) -> Callable:
    """Sets the duration and extremes for an event's pre- and post- trigger interval."""
    def _format_result(
        result: ChangepointMET,
        gti: GTI,
    ) -> Event:
        """Transforms a single focus results (times expressed as mets) into an event formatted like:
        (bkg_pres_start, bkg_pre_ends, event_starts, event_ends, bkg_post_start, bkg_post_end)
        """
        significance, changepoint, trigtime = result
        bkg_pre = compute_bkg_pre(
            trigtime,
            gti,
            pre_t,
            pre_delta,
        )
        bkg_post_start, bkg_post_end = compute_bkg_post(
            trigtime,
            gti,
            post_t,
            post_delta,
        )
        event_interval = changepoint, bkg_post_start
        # noinspection PyTypeChecker
        return Event(*bkg_pre, *event_interval, bkg_post_start, bkg_post_end)
    return _format_result


def format_results(
    results: dict[GTI, list[ChangepointMET]],
        intervals_duration_seconds: float,  # (binning / alpha)
        preinterval_ends_seconds: float,  # (binning * m)
        postinterval_start_seconds: float,  # binning * skip
) -> list[Event]:
    """Transforms the triggers changepoints to a list of events obeying the rule
    that an event must always be comprised in a GTI."""
    format_result = _format(
        intervals_duration_seconds,
        preinterval_ends_seconds,
        intervals_duration_seconds,
        postinterval_start_seconds,
    )
    return [
        format_result(changepoint, gti)
        for gti, cps in results.items()
        for changepoint in cps
    ]
