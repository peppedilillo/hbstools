from math import isclose
from pathlib import Path
from typing import Iterable

from _bisect import bisect_right
import numpy as np
import pandas as pd

from hbstools.io import read_event_files
from hbstools.io import read_gti_data
from hbstools.types import Dataset
from hbstools.types import Event
from hbstools.types import GTI
from hbstools.types import MET


def catalog(data_paths: Iterable[tuple[Path]]) -> Dataset:
    """Takes an unsorted collection of (data_path, gti_path) paths, and sorts
    it by the start time of the first GTI."""
    sorted_data_paths = sorted(data_paths, key=lambda d: read_gti_data(d[-1])[0].start)
    return [
        ((fp, gtip), gti)
        for fp, gtip in sorted_data_paths
        for gti in read_gti_data(gtip)
    ]


def _overlap(x: GTI, y: GTI, abs_tol: float) -> bool:
    """`x` overlaps `y` if `x` starts before or at the end of `y`"""
    assert x.start < y.end
    return isclose(x.end, y.start, abs_tol=abs_tol) or (y.start < x.end)


def _between(df, start_time: MET, end_time: MET) -> pd.DataFrame:
    return df[(df["TIME"] >= start_time) & (df["TIME"] < end_time)]


# this functions collates different datasets together, if the datasets are "adjacent"
# in time. the reason for doing this is to minimize the times in which we restart the
# trigger algorithm because at each restart a dead time is required to initialize an
# estimate of the background. other than the dataframe itself a gti is returned.
# this gtis annotates the data chunk's start and end time. this information is useful
# to perform further processing (e.g. formatting) of trigger events.
def stream(
    dataset: Dataset, abs_tol: float = 0.5
) -> Iterable[tuple[pd.DataFrame, GTI]]:
    """
    This function takes a dataset and returns an iterator, which will get you
    a (gti, DataFrame) tuple a time.
     The intended usage goes like:
    ```
    for gti, df in stream(dataset):
        run trigger on df

    :param dataset: a dataset, see `catalog`
    :param abs_tol: gtis closer than this are collated together
    :return:
    """
    # fp is for filepath, `p` prefix is for `pointer` or `past`
    ((pfp, _), pgti), *dataset = dataset
    df = read_event_files(pfp)
    pdf = _between(df, *pgti)
    for (fp, _), gti in dataset:
        # next line avoids multiple reads of the same data file.
        df = read_event_files(fp) if fp != pfp else df
        if not _overlap(pgti, gti, abs_tol):
            yield pdf, pgti
            pdf = _between(df, *gti)
            pgti = gti
        else:
            pdf = pd.concat((pdf, _between(df, max(pgti.end, gti.start), gti.end)))
            pgti = GTI(pgti.start, gti.end)
        pfp = fp
    yield pdf, pgti
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


def map_event_to_files(events: list[Event], dataset: Dataset) -> dict[Event, Path]:
    """Returns a dictionary of events and filepath. The filepath associated to
    an event is the one containing the event's start MET."""
    datafiles, gtis = list(zip(*dataset))
    gtis_starts = [gti.start for gti in gtis]
    ids = [bisect_right(gtis_starts, e.start) - 1 for e in events]
    assert -1 not in ids
    roots = [datafiles[i] for i in ids]
    return {event: root for event, root in zip(events, roots)}
