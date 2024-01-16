from pathlib import Path

from astropy.io import fits  # type: ignore[import-untyped]
import pandas as pd

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


def get_gtis(data_folder: str | Path) -> list[GTI]:
    def _get_gtis(gti_path: str | Path) -> list[GTI]:
        """Finds the GTIs."""
        gtis_hdul = fits.open(gti_path)
        return [GTI(start, stop) for start, stop in gtis_hdul[1].data]

    return _get_gtis(path_gtis(data_folder))


def get_merged_data(data_folder: str | Path) -> pd.DataFrame:
    def _get_merged_data(
        xdata_path: str | Path, sdata_path: str | Path
    ) -> pd.DataFrame:
        """Opens the data files and merges them"""
        xdata_df = pd.DataFrame(fits.open(xdata_path)[1].data)
        sdata_df = pd.DataFrame(fits.open(sdata_path)[1].data)
        return (
            pd.concat([xdata_df, sdata_df])
            .reset_index(drop=True)
            .sort_values(by=["TIME"])
        )

    return _get_merged_data(path_xdata(data_folder), path_sdata(data_folder))


def write_ttis_to_fits(ttis: pd.DataFrame, path: Path | str):
    """Write and save fits from tti dataframe"""
    content = fits.BinTableHDU.from_columns(ttis.to_records(index=False))
    content.writeto(path, overwrite=True)
