from math import isclose
from pathlib import Path
from typing import Iterator, Sequence

import pandas as pd

from hbstools.io import get_gtis
from hbstools.io import get_merged_data
from hbstools.types import GTI


def merge_overlapping_gtis(gtis: list[GTI], tolerance: float = 1.0) -> list[tuple]:
    """Merges CDI if they overlap. eg:
    F([(1, 2), (3, 4), (4, 5)]) = [(1, 2), (3, 5)]"""

    def overlap(x: GTI, y: GTI, abs_tol: float):
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
            yield from f(gtis, files[1:], concatenate(df, get_merged_data(new_file)))
        else:
            yield between(df, gti.start, gti.end), gti
            yield from f(gtis[1:], files, after(df, gti.end))

    all_gtis = [gti for gtis in [get_gtis(dp) for dp in data_folders] for gti in gtis]
    gtis = merge_overlapping_gtis(all_gtis)
    return f(gtis, data_folders[1:], get_merged_data(data_folders[0]))
