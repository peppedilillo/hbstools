from pathlib import Path
from typing import Iterator, Sequence

import numpy as np
import pandas as pd
from rich.console import Console
from rich.progress import track

from hbstools.data import get_data
from hbstools.io import get_gtis
from hbstools.poissonfocus import PoissonFocusDES
from hbstools.types import Change
from hbstools.types import ChangeMET
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


class Search:
    """Base search class."""

    def __init__(
        self,
        binning: float,
        skip: int,
        energy_lims: tuple[float, float],
        algorithm_params: dict,
        console: Console | None = None,
    ):
        self.binning = binning
        self.skip = skip
        self.enlims = energy_lims
        self.algorithm_params = algorithm_params
        self.algorithm = PoissonFocusDES
        self.console = console

    def __call__(self, dataset: Sequence[Path | str]) -> pd.DataFrame:
        return self.make_ttis(self.run_on_dataset(dataset))

    @staticmethod
    def make_bins(
        start: float,
        stop: float,
        step: float,
    ) -> np.ndarray:
        """Return bins including last one, which contains stop."""
        num_intervals = int((stop - start) / step + 1)
        return np.linspace(start, start + num_intervals * step, num_intervals + 1)

    def filter_data(
        self,
        data: pd.DataFrame,
    ):
        low, hi = self.enlims
        return data[(data["ENERGY"] >= low) & (data["ENERGY"] < hi)]

    def bin_data(
        self, data: pd.DataFrame, gti: GTI, binning: float
    ) -> tuple[np.ndarray, np.ndarray]:
        """Bins data in time and filters in energy."""
        bins = self.make_bins(gti.start, gti.end, binning)
        return np.histogram(data["TIME"], bins=bins)

    def run(
        self,
        xs: Sequence[int],
        bins: Sequence[float] | None = None,
    ) -> Change:
        """Run the algorithm one step a time."""
        algorithm = self.algorithm(**self.algorithm_params)
        for t, x_t in enumerate(xs):
            significance, offset = algorithm.step(x_t)
            if significance:
                return significance, t - offset + 1, t
        return 0.0, len(xs) + 1, len(xs)

    def run_on_segment(
        self,
        counts: np.ndarray,
        bins: np.ndarray,
    ) -> list[Change]:
        """Runs on binned data restarting the algorithm after sleep."""

        def f(cs, bs, skip, acc):
            if not len(cs):
                return []
            s, cp, tt = self.run(cs, bs)
            # if ended with no trigger algorithm returns 0, len(xs) + 1, len(xs)
            t = [(s, acc + cp, acc + tt)] if tt >= cp else []
            return t + f(cs[tt + skip :], bs[tt + skip :], skip, acc + tt + skip)

        return f(counts, bins, self.skip, 0)

    def run_on_dataset(
        self,
        dataset: Sequence[Path | str],
    ) -> dict[GTI, list[ChangeMET]]:
        """Runs on every gtis and returns anomalies."""

        def map_bins_to_times(rs, bs):
            """FOCuS time are expressed as bin-indeces.
            This transform from indeces to actual time."""
            return [(r[0], bs[r[1]], bs[r[2]]) for r in rs]

        def progressbar(gen: Iterator):
            return track(
                gen,
                description="[dim cyan]\[$Running]",
                transient=True,
                console=self.console,
            )

        results = {}
        sorted_folders = sorted(dataset, key=lambda d: get_gtis(d)[0].start)
        _get_data = (
            progressbar(get_data(sorted_folders))
            if self.console
            else get_data(sorted_folders)
        )
        for data, gti in _get_data:
            data = self.filter_data(data)
            counts, bins = self.bin_data(data, gti, self.binning)
            anomalies = self.run_on_segment(counts, bins)
            results[gti] = map_bins_to_times(anomalies, bins)

            # fmt: off
            if self.console:
                self.console.log(f"[dim]On GTI chunk {gti.start:.1f}-{gti.end:.1f}")
                if anomalies:
                    self.console.log(f"Found [b]{len(anomalies)}[/] transient{'s' if len(anomalies) > 1 else ''}")
                    for i, r in enumerate(results[gti], start=1):
                        self.console.log(f"[dim]MET time {r[2]:.1f}, GTI+{r[2] - gti.start:.1f}s-{r[2] - gti.start:.1f}s.[/]")
            # fmt: on
        return results

    def compute_bkg_pre(
        self,
        trigtime: MET,
        changepoint: MET,
        gti: GTI,
    ) -> tuple[float, float]:
        """Times from changepoint of an interval ending `m` steps behind the triggertime,
        with duration binning / alpha."""
        binning = self.binning
        duration = binning / self.algorithm_params["alpha"]
        m = self.algorithm_params["m"]

        if trigtime - m * binning - duration < gti.start:
            start = -(changepoint - gti.start)
            end = -m * binning + (trigtime - changepoint)
        else:
            end = (trigtime - changepoint) - binning * m
            start = end - duration
        return start, end

    def compute_bkg_post(
        self,
        trigtime: MET,
        changepoint: MET,
        gti: GTI,
    ) -> tuple[float, float]:
        """Times from changepoint of an interval starting after skip steps from triggertime,
        with duration binning / alpha."""
        binning = self.binning
        duration = binning / self.algorithm_params["alpha"]
        skip = self.skip

        if changepoint + (trigtime - changepoint) + binning * skip + duration > gti.end:
            end = gti.end - changepoint
            start = max(end - duration, binning)
        else:
            start = (trigtime - changepoint) + binning * skip
            end = start + duration
        return start, end

    def format_result(
        self,
        result: ChangeMET,
        gti: GTI,
    ) -> TTI:
        """Transforms focus results (times expressed as mets) into events formatted like:
        (bkg_pres_start, bkg_pre_ends, event_starts, event_ends, bkg_post_start, bkg_post_end)
        """
        significance, changepoint, trigtime = result
        bkg_pre = self.compute_bkg_pre(trigtime, changepoint, gti)
        bkg_post_start, bkg_post_end = self.compute_bkg_post(trigtime, changepoint, gti)
        event_interval = changepoint, changepoint + bkg_post_start
        return *bkg_pre, *event_interval, bkg_post_start, bkg_post_end

    def make_ttis(
        self,
        results: dict[GTI, list[ChangeMET]],
    ) -> pd.DataFrame:
        """Puts ttis into a dataframe"""
        formatted_results = []
        for gti in results.keys():
            for anomaly in results[gti]:
                formatted_results.append(self.format_result(anomaly, gti))
        return pd.DataFrame(formatted_results, columns=FCOLS)
