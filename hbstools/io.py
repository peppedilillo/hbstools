from pathlib import Path

from astropy.io import fits  # type: ignore[import-untyped]
import pandas as pd
from pandas.api.types import CategoricalDtype

from hbstools.types import GTI


def path_data(data_folder: str | Path) -> Path:
    """Returns path to x events"""
    return Path(data_folder).joinpath("out_lv1_cl.evt")


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
        data_path: str | Path
    ) -> pd.DataFrame:
        """Opens the data files and merges them"""
        with fits.open(data_path) as hdul:
            sdata_df = pd.DataFrame(hdul[1].data)
        category_quads_t = CategoricalDtype(categories=[0, 1, 2, 3], ordered=True)
        category_etype_t = CategoricalDtype(categories=[1, 2])
        return sdata_df.astype({
            "QUADID": category_quads_t,
            "EVTYPE": category_etype_t
        })

    return _read_event_files(path_data(data_folder))
