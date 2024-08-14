from pathlib import Path
from typing import Callable, Iterable

import pandas as pd
from rich.console import Console
from rich.progress import track

from hbstools.data import catalog
from hbstools.data import filter_energy
from hbstools.data import stream
from hbstools.format import format_results
import hbstools.trigger as trig
from hbstools.types import ChangepointMET
from hbstools.types import Dataset
from hbstools.types import GTI
from hbstools.types import Event


def search_set(
    binning: float,
    skip: int,
    energy_lims: tuple[float, float],
    algorithm_params: dict,
    log=lambda x: x,
):
    """Pepare a trigger algorithms and store parameters. Keep this pure."""

    def search_run(
        datastream: Iterable[tuple[pd.DataFrame, GTI]]
    ) -> dict[GTI, list[ChangepointMET]]:
        """Launch the search on every data table."""
        return {
            gti: run(filter_energy(data, energy_lims), gti) for data, gti in datastream
        }

    run = log(trig.trigger_df_set(binning, skip, algorithm_params))
    return search_run


def search_log(write: Callable):
    """This is where console writing, error checking and other housekeeping task happens.
    First, we take from the user a function which we will use for writing our messages.
    For example, we can take `console.log`, `print`, some function to write to file, or
    a combination of all."""

    def helper(f: Callable):
        """We take `f` the function which launches the algorithm on a single data table and
        return a troy horse. This troy horse looks exactly as `f` but intercepts its
        output and input, carrying out all the housekeeping tasks we need to perform."""

        def wrapper(data: pd.DataFrame, gti: GTI) -> list[ChangepointMET]:
            """The troy horse."""
            write(f"[dim]On GTI {gti.start:.0f}, {gti.end:.0f}..")
            try:
                rs = f(data, gti)
            except ValueError:
                write(f"[red]Error: invalid algorithm input.[/]")
                return []
            if rs:
                write(f"Found [b]{len(rs)}[/] transient{'s' if len(rs) > 1 else ''}")
                for r in rs:
                    write(f"[dim]|- MET {r[2]:1.1f}, (+{r[2] - gti.start:1.1f} s).[/]")
            return rs

        return wrapper

    return helper


def search(
    dataset: Dataset,
    configuration: dict,
    console: Console | None = None,
) -> list[Event]:
    """
    An interface to search.

    :param dataset: a dataset to be searched.
    :param configuration: an algorithm configuration.
    :param console: a rich console for writing.
    :return:
    """
    if console is not None:
        _log = search_log(console.log)
        datastream = track(
            stream(dataset),
            description="[dim cyan](Running..)",
            transient=True,
            console=console,
        )
    else:
        _log = search_log(lambda _: None)
        datastream = stream(dataset)

    algorithm_params = configuration["algorithm_params"]
    search_run = search_set(
        configuration["binning"],
        configuration["skip"],
        configuration["energy_lims"],
        configuration["algorithm_params"],
        _log,
    )
    results = search_run(datastream)
    return format_results(
        results=results,
        intervals_duration_seconds=configuration["binning"] / algorithm_params["alpha"],
        preinterval_ends_seconds=configuration["binning"] * algorithm_params["m"],
        postinterval_start_seconds=configuration["binning"] * configuration["skip"],
    )
