from pathlib import Path
from typing import Callable, Iterable

import pandas as pd
from rich.console import Console
from rich.progress import track

import hbstools.trigger as trig
from hbstools.data import catalog, filter_energy, stream
from hbstools.format import format_results
from hbstools.types import GTI, ChangepointMET


def search_set(
    binning: float,
    skip: int,
    energy_lims: tuple[float, float],
    algorithm_params: dict,
    log=lambda x: x,
):
    """Pepare a trigger algorithms and store parameters. Keep this pure."""

    def search_run(
        data_stream: Iterable[tuple[pd.DataFrame, GTI]]
    ) -> dict[GTI, list[ChangepointMET]]:
        """Launch the search on every data table."""
        return {
            gti: run(filter_energy(data, energy_lims), gti) for data, gti in data_stream
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
            (
                write(f"Found [b]{len(rs)}[/] transient{'s' if len(rs) > 1 else ''}")
                if rs
                else None
            )
            for r in rs:
                write(f"[dim]|- MET {r[2]:1.1f}, (+{r[2] - gti.start:1.1f} s).[/]")
            return rs

        return wrapper

    return helper


def search(
    data_folders: Iterable[Path | str],
    configuration: dict,
    console: Console | None = None,
):
    """
    An interface to search.

    :param data_folders: an iterable for the inputs. file ordering is not assumed.
    :param configuration: an algorithm configuration.
    :param console: a rich console for writing.
    TODO: change console to a Callable.
    :return:
    """
    dataset = catalog(data_folders)
    if console is not None:
        _log = search_log(console.log)
        data_stream = track(
            stream(dataset),
            description="[dim cyan](Running..)",
            transient=True,
            console=console,
        )
    else:
        _log = search_log(lambda _: None)
        data_stream = stream(dataset)

    algorithm_params = configuration["algorithm_params"]
    return format_results(
        results=search_set(
            configuration["binning"],
            configuration["skip"],
            configuration["energy_lims"],
            configuration["algorithm_params"],
            _log,
        )(
            data_stream,
        ),
        intervals_duration_seconds=configuration["binning"] / algorithm_params["alpha"],
        preinterva_ends_seconds=configuration["binning"] * algorithm_params["m"],
        postinterval_start_seconds=configuration["binning"] * configuration["skip"],
    )
