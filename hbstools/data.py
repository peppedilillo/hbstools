from functools import reduce
from itertools import pairwise
from math import isclose
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from hbstools.io import read_event_files, read_gti_files
from hbstools.types import GTI, MET, Dataset


def catalog(data_folders: Iterable[Path | str]) -> Dataset:
    """Takes an unordered collection of paths, orders them and returns a Dataset."""
    unsorted_gtis = {dp: read_gti_files(dp) for dp in data_folders}
    sorted_folders = sorted(data_folders, key=lambda d: unsorted_gtis[d][0].start)
    sorted_gtis = [unsorted_gtis[dp] for dp in sorted_folders]

    def is_sorted(xs):
        return len(xs) < 2 or reduce(lambda acc, x: acc and x[0] < x[1], pairwise(xs))

    assert is_sorted([gti.start for gtis in sorted_gtis for gti in gtis]) and is_sorted(
        [gti.end for gtis in sorted_gtis for gti in gtis]
    )
    return [(gti, dp) for gtis, dp in zip(sorted_gtis, sorted_folders) for gti in gtis]


# this functions collates different datasets together, if the datasets are "adjacent"
# in time. the reason for doing this is to minimize the times in which we restart the
# trigger algorithm because at each restart a dead time is required to initialize an
# estimate of the background. other than the dataframe itself a gti is returned.
# this gtis annotates the data chunk's start and end time. this information is useful
# to perform further processing (e.g. formatting) of trigger events.
def stream(
    dataset: Dataset, abs_tol: float = 0.5
) -> Iterable[tuple[GTI, pd.DataFrame]]:
    """This function takes a dataset and returns an iterator, which will get you
    a (gti, DataFrame) tuple a time. The intended usage goes like:
    ```
    for gti, df in stream(dataset):
        run trigger on df
    ```
    """

    def overlap(x: GTI, y: GTI, abs_tol: float) -> bool:
        """`x` overlaps `y` if `x` starts before or at the end of `y`"""
        assert x.start < y.end
        return isclose(x.end, y.start, abs_tol=abs_tol) or (y.start < x.end)

    def between(df, start_time: MET, end_time: MET) -> pd.DataFrame:
        return df[(df["TIME"] >= start_time) & (df["TIME"] < end_time)]

    (last_gti, data_folder), *dataset = dataset
    last_df = between(read_event_files(data_folder), *last_gti)
    for gti, data_folder in dataset:
        df = read_event_files(data_folder)
        if not overlap(last_gti, gti, abs_tol):
            yield last_df, last_gti
            last_df = between(df, *gti)
            last_gti = gti
        else:
            last_df = pd.concat(
                (last_df, between(df, max(last_gti.end, gti.start), gti.end))
            )
            last_gti = GTI(last_gti.start, gti.end)
    yield last_df, last_gti
    return


def filter_energy(
    data: pd.DataFrame,
    energy_lims: tuple[float, float],
) -> pd.DataFrame:
    """Filters data in an energy band."""
    low, hi = energy_lims
    return data[(data["ENERGY"] >= low) & (data["ENERGY"] < hi)]


def _histogram(
    data: pd.Series,
    start: float,
    stop: float,
    binning: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Bins data in time and returns counts and bins. The `counts` array has
    length int((stop - start) / binning + 1), the `bins` array is
    longer by one unit. The last element of `bins` is guaranteed to be
    greater-equal than stop."""
    num_intervals = int((stop - start) / binning + 1)
    counts, bins = np.histogram(
        data,
        range=(start, start + num_intervals * binning),
        bins=num_intervals,
    )
    return counts, bins


def histogram(
    data: pd.DataFrame,
    gti: GTI,
    binning: float,
) -> tuple[np.ndarray, np.ndarray]:
    """A specialized histogram to histogram time event lists."""
    return _histogram(data["TIME"], *gti, binning)


def histogram_quadrants(
    data: pd.DataFrame,
    gti: GTI,
    binning: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Bins data in time, separating data from different quadrants."""
    _, bins = _histogram(pd.Series([]), *gti, binning)
    quadrant_counts = [
        histogram(quadrant_data, gti, binning)[0]
        for _, quadrant_data in data.groupby("QUADID", observed=False)
    ]
    return np.vstack(quadrant_counts), bins
