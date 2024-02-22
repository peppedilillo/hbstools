from pathlib import Path

from astropy.io import fits  # type: ignore[import-untyped]
import pandas as pd
from pandas.api.types import CategoricalDtype

from hbstools.types import GTI


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
        """Finds the GTIs."""
        gtis_hdul = fits.open(gti_path)
        return [GTI(start, stop) for start, stop in gtis_hdul[1].data]

    return _read_gti_files(path_gtis(data_folder))


def read_event_files(data_folder: str | Path) -> pd.DataFrame:
    def _read_event_files(
        xdata_path: str | Path, sdata_path: str | Path
    ) -> pd.DataFrame:
        """Opens the data files and merges them"""
        xdata_df = pd.DataFrame(fits.open(xdata_path)[1].data)
        sdata_df = pd.DataFrame(fits.open(sdata_path)[1].data)
        category_quads_t = CategoricalDtype(categories=[0,1,2,3], ordered=True)
        return (
            pd.concat([xdata_df, sdata_df])
            .reset_index(drop=True)
            .sort_values(by=["TIME"])
            .astype({'QUADID': category_quads_t})
        )

    return _read_event_files(path_xdata(data_folder), path_sdata(data_folder))


def write_ttis_to_fits(ttis: pd.DataFrame, path: Path | str):
    """Write and save fits from tti dataframe"""
    content = fits.BinTableHDU.from_columns(ttis.to_records(index=False))
    content.writeto(path, overwrite=True)
