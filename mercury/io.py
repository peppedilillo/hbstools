from datetime import datetime
from math import log10
from pathlib import Path
from uuid import uuid4

from astropy.io import fits
import numpy as np
import pandas as pd
import yaml

import hbstools as hbs
from hbstools.data import map_event_to_files
from hbstools.types import Dataset
from hbstools.types import Event


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
    index = {"uuid": uuid4().hex, "mappings": (fmap := {})}
    pad = int(log10(len(events))) + 1  # for filename padding
    for n, (event, (_, gti_path)) in enumerate(
        map_event_to_files(events, dataset).items()
    ):
        gti_content = hbs.io.read_gti_file(gti_path)
        _write_src(
            event,
            src_path := dir_path / f"event-src-{n:0{pad}}.fits",
            configuration,
            gti_content,
        )
        _write_bkg(
            event,
            bkg_path := dir_path / f"event-bkg-{n:0{pad}}.fits",
            configuration,
            gti_content,
        )
        fmap[src_path.name] = {"root": str(gti_path.parent.absolute()), "type": "src"}
        fmap[bkg_path.name] = {"root": str(gti_path.parent.absolute()), "type": "bkg"}

    with open(dir_path / index_fname, "w") as f:
        yaml.dump(index, f)


def _write_src(
    event: Event,
    filepath: Path,
    configuration: dict,
    gti: tuple[np.recarray, fits.Header] | None = None,
):
    """A helper for writing an event's source output file.
    Primary HDU header is the HDU of the GTI where the event start time is."""
    data = pd.DataFrame(
        {
            "START": [
                event.start,
            ],
            "STOP": [
                event.stop,
            ],
        }
    ).to_records(index=False)
    primary = fits.PrimaryHDU(
        header=fits.Header(_tag({}).items()),
    )
    gti_data, gti_header = gti
    gti = fits.BinTableHDU.from_columns(
        gti_data,
        header=gti_header,
        name="STDGTI",
    )
    data = fits.BinTableHDU.from_columns(
        data,
        header=fits.Header(_compile_data_header(configuration).items()),
        name="TRIGGERS",
    )
    fits.HDUList([primary, gti, data]).writeto(filepath)


def _write_bkg(
    event: Event,
    filepath: Path,
    configuration: dict,
    gti: tuple[np.recarray, fits.Header] | None = None,
):
    """A helper for writing an event's background output file.
    Primary HDU header is the HDU of the GTI where the event start time is."""
    data = pd.DataFrame(
        {
            "BKG_START": [event.bkg_pre_start, event.bkg_post_start],
            "BKG_STOP": [event.bkg_post_start, event.bkg_post_stop],
        }
    ).to_records(index=False)
    primary = fits.PrimaryHDU(
        header=fits.Header(_tag({}).items()),
    )
    gti_data, gti_header = gti
    gti = fits.BinTableHDU.from_columns(
        gti_data,
        header=gti_header,
        name="STDGTI",
    )
    data = fits.BinTableHDU.from_columns(
        data,
        header=fits.Header(_compile_data_header(configuration).items()),
        name="TRIGGERS",
    )
    fits.HDUList([primary, gti, data]).writeto(filepath)


def write_catalog(
    events: list[Event],
    filepath: Path | str,
    configuration: dict,
):
    """Write results to fits under catalog mode."""
    data = pd.DataFrame(
        {
            "BKG_PRE_START": [e.bkg_pre_start for e in events],
            "BKG_PRE_STOP": [e.bkg_pre_stop for e in events],
            "START": [e.start for e in events],
            "STOP": [e.stop for e in events],
            "BKG_POST_START": [e.bkg_post_start for e in events],
            "BKG_POST_STOP": [e.bkg_post_stop for e in events],
        }
    ).to_records(index=False)
    primary = fits.PrimaryHDU(
        header=fits.Header(_tag({}).items()),
    )
    data = fits.BinTableHDU.from_columns(
        data,
        header=fits.Header(_compile_data_header(configuration).items()),
        name="TRIGGERS",
    )
    fits.HDUList([primary, data]).writeto(filepath)


def _flat(d: dict) -> dict:
    """Flattens a dictionary recursively.
    Example:
    In  = {"a": 1, "b": 2, "c": {"d": 3, "e": 4}}
    Out = {"a": 1, "b": 2, "d": 3, "e": 4}
    """

    def helper(s, acc):
        if not s:
            return acc
        k, v = s.popitem()
        if not isinstance(v, dict):
            return helper(s, acc | {k: v})
        else:
            return helper(s, helper(v, {}))

    return helper(d, {})


def _short(d: dict) -> dict:
    """Avoids the creation of an un-standard header card"""
    return {k[:8]: v for k, v in d.items()}


def _vstring(d: dict) -> dict:
    """Writing strings we can ignore iterable values in configurations."""
    return {k: str(v) for k, v in d.items()}


def _tag(d: dict) -> dict:
    """Adds a number of mercury-specific tags"""
    # TODO: add a version tag.
    return d | {
        "CREATOR": "hbst-mercury",
        "DATE": datetime.now().strftime("%m/%d/%Y-%H:%M:%S"),
    }


def _compile_data_header(configuration: dict) -> dict:
    """Starting from a configuration returns a dictionary which can be used
    as a FITS header."""
    return _vstring(_short(_tag(_flat(configuration))))
