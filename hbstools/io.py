from functools import cache
from pathlib import Path

import numpy as np
from astropy.io import fits  # type: ignore[import-untyped]
from astropy.table import Table
import pandas as pd
from pandas.api.types import CategoricalDtype

from hbstools.types import GTI


# GTI files are small and we need to open them a few times
@cache
def read_gti_file(gti_path: str | Path) -> tuple[np.recarray, fits.Header]:
    """Reads a GTI and caches its content."""
    return fits.getdata(gti_path, header=True)


def read_gti_data(gti_path: str | Path) -> list[GTI]:
    """Returns a list of GTI from file path."""
    data, _ = read_gti_file(gti_path)
    return [GTI(*map(float, content)) for content in data]


def read_gti_header(gti_path: str | Path) -> fits.Header:
    """Returns a GTI header."""
    _, header = read_gti_file(gti_path)
    return header


def read_event_files(data_path: str | Path) -> pd.DataFrame:
    """Opens the data files and merges them"""
    data_df = Table.read(data_path, hdu=1, format="fits").to_pandas()
    category_quads_t = CategoricalDtype(categories=[0, 1, 2, 3], ordered=True)
    return data_df.astype(
        {
            "QUADID": category_quads_t,
        }
    )
