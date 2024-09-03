from datetime import datetime
from math import log10
from pathlib import Path
from uuid import uuid4

import numpy as np
import pandas as pd
import yaml
from astropy.io import fits

import hbstools as hbs
from hbstools.data import map_event_to_files
from hbstools.types import Event, Dataset


def write_library(
        events: list[Event],
        dataset: Dataset,
        configuration: dict,
        dir_path: Path,
        index_fname: str,
):
    """
    Writes events to multiple fits files and store an index of them which
    can be used to merge the results to the dataset.

    :param events:
    :param dataset:
    :param configuration:
    :param dir_path: shall point to a directory. both the events and the index
    are given default names.
    :param index_fname: name of the index file
    """

    def write_src(event: Event, filepath: Path, header: fits.Header | None = None):
        """A helper for writing an event's source output file.
        Primary HDU header is the HDU of the GTI where the event start time is."""
        data = pd.DataFrame([event])[["start", "end"]].to_records(index=False)
        _write_results_fits(
            data,
            filepath,
            primary_header=header,
            data_header=fits.Header(_compile_data_header(configuration).items()),
        )

    def write_bkg(e: Event, filepath: Path, header: fits.Header | None = None):
        """A helper for writing an event's background output file.
        Primary HDU header is the HDU of the GTI where the event start time is."""
        data = pd.DataFrame(
            {
                "bkg_start": [e.bkg_pre_start, e.bkg_post_start],
                "bkg_end": [e.bkg_post_start, e.bkg_post_end],
            }
        ).to_records(index=False)
        _write_results_fits(
            data,
            filepath,
            primary_header=header,
            data_header=fits.Header(_compile_data_header(configuration).items()),
        )

    index = {"uuid": uuid4().hex, "mappings": (fmap := {})}
    pad = int(log10(len(events))) + 1  # for filename padding
    for n, (event, (_, gti_path)) in enumerate(
        map_event_to_files(events, dataset).items()
    ):
        header = hbs.io.read_gti_header(gti_path)
        write_src(event, src_path := dir_path / f"event-src-{n:0{pad}}.fits", header)
        write_bkg(event, bkg_path := dir_path / f"event-bkg-{n:0{pad}}.fits", header)
        fmap[src_path.name] = {"root": str(gti_path.parent.absolute()), "type": "src"}
        fmap[bkg_path.name] = {"root": str(gti_path.parent.absolute()), "type": "bkg"}

    with open(dir_path / index_fname, "w") as f:
        yaml.dump(index, f)


def write_catalog(
        events: list[Event],
        filepath: Path | str,
        configuration: dict,
):
    """Write results to fits under catalog mode."""
    data = pd.DataFrame(events).to_records(index=False)
    _write_results_fits(
        data,
        filepath,
        primary_header=fits.Header(_tag({}).items()),
        data_header=fits.Header(_compile_data_header(configuration).items()),
    )


def _write_results_fits(
        data: np.recarray,
        filepath: Path | str,
        primary_header: fits.Header | None = None,
        data_header: fits.Header | None = None,
):
    """Savesto a single fits file writing configuration in the file data HDU header."""
    primary = fits.PrimaryHDU(
        header=primary_header,
    )
    data = fits.BinTableHDU.from_columns(
        data,
        header=data_header,
        name="TRIGGERS",
    )
    fits.HDUList([primary, data]).writeto(filepath)


def _flat(d: dict):
    """Flattens a dictionary recursively.
    Example:
    In  = {"a": 1, "b": 2, "c": {"d": 3, "e": 4}}
    Out = {"a": 1, "b": 2, "d": 3, "e": 4}
    """
    def helper(d, acc):
        if not d:
            return acc
        k, v = d.popitem()
        if not isinstance(v, dict):
            return helper(d, acc | {k: v})
        else:
            return helper(d, helper(v, {}))

    return helper(d, {})


def _short(d: dict):
    """Avoids the creation of an un-standard header card"""
    return {k[:8]: v for k, v in d.items()}


def _vstring(d: dict):
    """Writing strings we can ignore iterable values in configurations."""
    return {k: str(v) for k, v in d.items()}


def _tag(d: dict):
    """Adds a number of mercury-specific tags"""
    # TODO: add a version tag.
    return d | {
        "CREATOR": "hbst-mercury",
        "DATE": datetime.now().strftime("%m/%d/%Y-%H:%M:%S"),
    }


def _compile_data_header(configuration: dict):
    return _vstring(_short(_tag(_flat(configuration))))
