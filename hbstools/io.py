from pathlib import Path

import pandas as pd
from astropy.io import fits  # type: ignore[import-untyped]
from pandas.api.types import CategoricalDtype

from hbstools.types import GTI, Event


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
