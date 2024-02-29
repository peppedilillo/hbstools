from pathlib import Path
from typing import Iterator, Sequence

import pandas as pd
from rich.console import Console
from rich.progress import track

from hbstools.data import filter_energy
from hbstools.data import get_data
import hbstools.trigger as trig
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
        self.energy_lims = energy_lims
        self.algorithm_params = algorithm_params
        self.algorithm_type = trig.match_algorithm(algorithm_params)
        self.console = console

    def __call__(self, dataset: Sequence[Path | str]) -> pd.DataFrame:
        return self.make_ttis(self.run_on_dataset(dataset))

    def run_on_dataset(
        self,
        folders: Sequence[Path | str],
    ) -> dict[GTI, list[ChangepointMET]]:
        """Runs on every gtis and returns anomalies."""
        def pbar(gen: Iterator):
            """A progress bar."""
            return track(
                gen,
                description="[dim cyan](Running..)",
                transient=True,
                console=self.console,
            )

        _get_data = pbar(get_data(folders)) if self.console else get_data(folders)
        results = {}
        run = trig.set(self.algorithm_params)
        for data, gti in _get_data:
            if self.console is not None:
                self.console.log(f"[dim]On GTI {gti.start:.0f}, {gti.end:.0f}..")
            filtered_data = filter_energy(data, self.energy_lims)
            try:
                anomalies = run(filtered_data, gti, self.binning, self.skip)
            except ValueError:
                if self.console is not None:
                    self.console.log(f"[red]Error: invalid algorithm input.[/]")
                continue

            results[gti] = anomalies

            # fmt: off
            if self.console is not None and anomalies:
                self.console.log(
                    f"Found [b]{len(anomalies)}[/] transient{'s' if len(anomalies) > 1 else ''}")
                for i, r in enumerate(results[gti], start=1):
                    self.console.log(
                        f"[dim]|- MET {r[2]:1.1f}, "
                        f"(+{r[2] - gti.start:1.1f} s).[/]"
                    )
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
        result: ChangepointMET,
        gti: GTI,
    ) -> TTI:
        """Transforms focus results (times expressed as mets) into events formatted like:
        (bkg_pres_start, bkg_pre_ends, event_starts, event_ends, bkg_post_start, bkg_post_end)
        """
        significance, changepoint, trigtime = result
        bkg_pre = self.compute_bkg_pre(trigtime, changepoint, gti)
        bkg_post_start, bkg_post_end = self.compute_bkg_post(trigtime, changepoint, gti)
        event_interval = changepoint, changepoint + bkg_post_start
        # noinspection PyTypeChecker
        return *bkg_pre, *event_interval, bkg_post_start, bkg_post_end

    def make_ttis(
        self,
        results: dict[GTI, list[ChangepointMET]],
    ) -> pd.DataFrame:
        """Puts ttis into a dataframe"""
        formatted_results = []
        for gti in results.keys():
            for anomaly in results[gti]:
                formatted_results.append(self.format_result(anomaly, gti))
        return pd.DataFrame(formatted_results, columns=FCOLS)
