from _bisect import bisect_right
from pathlib import Path

import pandas as pd
from astropy.io import fits  # type: ignore[import-untyped]
from pandas.api.types import CategoricalDtype
import yaml

from hbstools.types import GTI, Event, Dataset



def path_xdata(data_folder: str | Path) -> Path:
    """Returns path to x events"""
    return Path(data_folder).joinpath("out_x_cl.evt")


def path_sdata(data_folder: str | Path) -> Path:
    """Returns path to s events"""
    return Path(data_folder).joinpath("out_s_cl.evt")


def path_gtis(data_folder: str | Path) -> Path:
    """Returns path to gtis"""
    return Path(data_folder).joinpath("gti.fits")


def read_gti_files(data_folder: str | Path) -> list[GTI]:
    def _read_gti_files(gti_path: str | Path) -> list[GTI]:
        """Returns a list of GTI from file path."""
        with fits.open(gti_path) as hdul:
            return [GTI(*map(float, content)) for content in hdul[1].data]

    return _read_gti_files(path_gtis(data_folder))


def read_event_files(data_folder: str | Path) -> pd.DataFrame:
    def _read_event_files(
        xdata_path: str | Path, sdata_path: str | Path
    ) -> pd.DataFrame:
        """Opens the data files and merges them"""
        with fits.open(xdata_path) as hdul:
            xdata_df = pd.DataFrame(hdul[1].data)
        with fits.open(sdata_path) as hdul:
            sdata_df = pd.DataFrame(hdul[1].data)
        category_quads_t = CategoricalDtype(categories=[0, 1, 2, 3], ordered=True)
        return (
            pd.concat([xdata_df, sdata_df])
            .sort_values(by=["TIME"])
            .astype({"QUADID": category_quads_t})
            .reset_index(drop=True)
        )

    return _read_event_files(path_xdata(data_folder), path_sdata(data_folder))


def write_catalog(events: list[Event], filepath: Path | str, header: fits.Header | None = None):
    """Save an event list to a single fits file"""
    content = pd.DataFrame(events).to_records(index=False)
    fits.writeto(filename=filepath, data=content, header=header, overwrite=True)


def write_source_fits(event: Event, filepath: Path, header: fits.Header | None = None):
    content = pd.DataFrame([event])[["start", "end"]].to_records(index=False)
    fits.writeto(filename=filepath, data=content, header=header, overwrite=True)


def write_bkg_fits(event: Event, filepath: Path, header: fits.Header | None = None):
    content = pd.DataFrame(
        {
            "bkg_start": [event.bkg_pre_start, event.bkg_post_start],
            "bkg_end": [event.bkg_post_start, event.bkg_post_end],
        }
    ).to_records(index=False)
    fits.writeto(filename=filepath, data=content, header=header, overwrite=True)


def map_event_to_files(events: list[Event], dataset: Dataset) -> dict[Event, Path]:
    """Returns a dictionary of events and filepath. The filepath associated to
    an event is the one containing the event's start MET."""
    gtis, datafiles = list(zip(*dataset))
    gtis_starts = [gti.start for gti in gtis]
    ids = [bisect_right(gtis_starts, e.start) - 1 for e in events]
    assert -1 not in ids
    roots = [datafiles[i] for i in ids]
    return {event: root for event, root in zip(events, roots)}


INDEX_FILENAME = ".mercury-index.yaml"


def write_index(index: dict, path: Path):
    with open(path / INDEX_FILENAME, 'w') as f:
        yaml.dump(index, f)


def write_library(events: list[Event], dataset: Dataset, path: Path):
    from math import log10

    event_map = map_event_to_files(events, dataset)
    index = {}
    width = int(log10(len(events))) + 1  # filename padding
    for num, (event, root) in enumerate(event_map.items()):
        _, header = fits.getdata(path_gtis(root), header=True)
        suffix = f"-{num:0{width}}"

        write_source_fits(event, src_path := path / f"event-src{suffix}.fits", header)
        write_bkg_fits(event, bkg_path := path / f"event-bkg{suffix}.fits", header)
        index[src_path.name] = {
            "root": str(root.absolute()),
            "type": "src",
        }
        index[bkg_path.name] = {
            "root": str(root.absolute()),
            "type": "bkg",
        }

    write_index(index, path)
