"""
Ctypes binding to the header of the C implementation of Poisson FOCuS
with single exponential smoothing background estimate.
"""

import ctypes
import enum
from typing import Callable

import numpy as np
import numpy.typing as npt

from hbstools.triggers import _LIBCFOCUS
from hbstools.types import Changepoint

clib_pfs = ctypes.CDLL(_LIBCFOCUS)


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


class PoissonFocusSesCwrapper:
    """A wrapper to the C implementation of FOCuS with simple exp. smoothing."""

    def __init__(
        self,
        threshold_std: float,
        mu_min: float,
        alpha: float,
        m: int,
        sleep: int,
    ):
        self.check_init_parameters(threshold_std, mu_min, alpha, m, sleep)
        self.threshold_std = threshold_std
        self.mu_min = mu_min
        self.alpha = alpha
        self.m = m
        self.sleep = sleep
        self._call = self.bind_pfs_interface()

    def __call__(self, xs: npt.NDArray) -> Changepoint:
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
            self.sleep,
        )
        match error_code:
            case _Errors.INVALID_INPUT:
                raise ValueError("The inputs contain invalid entries.")
            case _Errors.INVALID_ALLOCATION:
                raise BufferError("Can't allocate memory for the algorithm.")
        return c.significance_std, c.changepoint, c.triggertime

    @staticmethod
    def check_init_parameters(
        threshold_std: float,
        mu_min: float,
        alpha: float,
        m: int,
        sleep: int,
    ):
        """Checks validity of initialization arguments."""
        check = PoissonFocusSesCwrapper.bind_pfs_check_init_parameters()
        error_code = check(
            threshold_std,
            mu_min,
            alpha,
            m,
            sleep,
        )
        match error_code:
            case _Errors.INVALID_INPUT:
                raise ValueError("The inputs contain invalid entries.")
        return

    @staticmethod
    def bind_pfs_check_init_parameters() -> Callable:
        pfs_check_init_parameters = clib_pfs.pfs_check_init_parameters
        pfs_check_init_parameters.restype = _Errors
        pfs_check_init_parameters.argtypes = [
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int,
            ctypes.c_int,
        ]
        return pfs_check_init_parameters

    @staticmethod
    def bind_pfs_interface() -> Callable:
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
            ctypes.c_int,
        ]
        return pfs_interface
