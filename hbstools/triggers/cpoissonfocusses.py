import ctypes
import enum
from typing import Callable

import numpy as np
import numpy.typing as npt


def library_path(libname: str) -> str:
    """Gets the path to the library, findind the appropriate extension based on
     OS. The lib is installed in `lib/pythonX.XX/site-packages/.sharedlibs`.
     Assumes this script to be installed at `lib/pythonX.XX/site-packages`"""
    from pathlib import Path
    import sys

    if sys.platform.startswith('win32'):
        suffix = ".dll"
    elif sys.platform.startswith('linux'):
        suffix = ".so"
    elif sys.platform.startswith('darwin'):
        suffix = ".dylib"
    else:
        raise OSError("System not supported")

    return str(Path(__file__).parent / f"{libname}{suffix}")


clib_pfs = ctypes.CDLL(library_path(".sharedlibs/lib-pfocus"))


class _Errors(enum.IntEnum):
    """Errors returned by the C implementation."""
    NO_ERRORS = 0
    INVALID_ALLOCATION = 1
    INVALID_INPUT = 2


class _Changepoint(ctypes.Structure):
    """This is a wrapper to our C definition of a changepoint."""
    _fields_ = [
        ("significance_std", ctypes.c_double),
        ("changepoint", ctypes.c_size_t),
        ("triggertime", ctypes.c_size_t),
    ]


class PoissonFocusSES:
    """A wrapper to the C implementation of FOCuS with simple exp. smoothing."""
    def __init__(
            self,
            threshold_std: float,
            mu_min: float,
            alpha: float,
            beta,  # TODO: remove this when adding dispatcher
            t_max,  # TODO: same as above
            m: int,
            sleep: int
    ):
        self.threshold_std = threshold_std
        self.mu_min = mu_min
        self.alpha = alpha
        self.m = m
        self.sleep = sleep
        self._call = self.bind()

    def run(self, xs: npt.NDArray[np.int_]) -> tuple[float, int, int]:
        c = _Changepoint()
        xs_length = len(xs)
        xs_pointer = xs.ctypes.data_as(ctypes.POINTER(ctypes.c_long))
        error_code = self._call(
            ctypes.byref(c),
            xs_pointer,
            xs_length,
            self.threshold_std,
            self.mu_min,
            self.alpha,
            self.m,
            self.sleep
        )

        match error_code:
            case _Errors.INVALID_INPUT:
                raise ValueError("The inputs contain invalid entries.")
            case _Errors.INVALID_ALLOCATION:
                raise BufferError("Can't allocate memory for the algorithm.")
        return c.significance_std, c.changepoint, c.triggertime

    @staticmethod
    def bind() -> Callable:
        pfs_interface = clib_pfs.pfs_interface
        pfs_interface.restype = _Errors
        pfs_interface.argtypes = [
            ctypes.POINTER(_Changepoint),
            ctypes.POINTER(ctypes.c_long),
            ctypes.c_size_t,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int,
            ctypes.c_int
        ]
        return pfs_interface
