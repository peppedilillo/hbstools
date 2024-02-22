from math import isclose
from pathlib import Path
from typing import Iterator, Sequence

import pandas as pd
import numpy as np

from hbstools.io import read_gti_files
from hbstools.io import read_event_files
from hbstools.types import GTI


def merge_overlapping_gtis(gtis: list[GTI], tolerance: float = 1.0) -> list[tuple]:
    """Merges CDI if they overlap. eg:
    F([(1, 2), (3, 4), (4, 5)]) = [(1, 2), (3, 5)]"""

    def overlap(x: GTI, y: GTI, abs_tol: float):
        """`x` overlaps `y` if `x` starts before or at the end of `y` """
        assert x.start < y.end
        return isclose(x.end, y.start, abs_tol=abs_tol) or (y.start < x.end)

    def f(xs: Sequence[GTI]):
        if len(xs) == 2:
            first, second = xs
            if overlap(first, second, tolerance):
                return [GTI(first.start, second.end)]
        if len(xs) > 2:
            first, second, *cdr = xs
            if overlap(first, second, tolerance):
                return f([GTI(first.start, second.end)] + cdr)
            else:
                return [first] + f([second, *cdr])
        return xs

    return f(gtis)


def get_data(data_folders: Sequence[Path | str]) -> Iterator[tuple[pd.DataFrame, GTI]]:
    """A generator which will get you one GTI dataframe a time."""

    def after(df, time):
        return df[df["TIME"] >= time]

    def before(df, time):
        return df[df["TIME"] < time]

    def between(df, start_time, end_time):
        return before(after(df, start_time), end_time)

    def concatenate(df1, df2):
        first_time = df2["TIME"].iloc[0]
        return pd.concat((before(df1, first_time), df2))

    def f(gtis, files, df):
        if len(gtis) == 0:
            return
        gti, *_ = gtis
        if gti.end > df["TIME"].iloc[-1]:
            if not files:
                yield after(df, gti.start), gti
                return
            new_file, *_ = files
            yield from f(gtis, files[1:], concatenate(df, read_event_files(new_file)))
        else:
            yield between(df, gti.start, gti.end), gti
            yield from f(gtis[1:], files, after(df, gti.end))

    all_gtis = [gti for gtis in [read_gti_files(dp) for dp in data_folders] for gti in gtis]
    gtis = merge_overlapping_gtis(all_gtis)
    return f(gtis, data_folders[1:], read_event_files(data_folders[0]))


def make_bins(
    start: float,
    stop: float,
    step: float,
) -> np.ndarray:
    """Return bins including last one, which contains stop."""
    num_intervals = int((stop - start) / step + 1)
    return np.linspace(start, start + num_intervals * step, num_intervals + 1)


def filter_energy(
    energy_lims: tuple[float, float],
    data: pd.DataFrame,
):
    """Filters data in an energy band."""
    low, hi = energy_lims
    return data[(data["ENERGY"] >= low) & (data["ENERGY"] < hi)]


def histogram(
    data: pd.DataFrame,
    gti: GTI,
    binning: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Bins data in time."""
    bins = make_bins(gti.start, gti.end, binning)
    counts, _ = np.histogram(data["TIME"], bins=bins)
    return counts, bins


def histogram_quadrants(
    data: pd.DataFrame,
    gti: GTI,
    binning: float,
) -> tuple[dict, np.ndarray]:
    """Bins data in time, separating data from different quadrants."""
    bins = make_bins(gti.start, gti.end, binning)
    return {
        quadid: histogram(quadrant_data, gti, binning)[0]
        for quadid, quadrant_data in data.groupby("QUADID")
    }, bins
