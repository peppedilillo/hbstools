"""
Ctypes binding to the header of the BFT C implementation.
"""

import ctypes
import enum
from typing import Callable

import numpy as np
import numpy.typing as npt

from hbstools.triggers import _LIBCFOCUS
from hbstools.triggers import TriggerAlgorithm
from hbstools.types import Changepoint

clib_bft = ctypes.CDLL(_LIBCFOCUS)

NUMDETECTOR = 4


class _Errors(enum.IntEnum):
    """Errors returned by the C implementation."""

    NO_ERRORS = 0
    INVALID_ALLOCATION = 1
    INVALID_INPUT = 2


class _Changepoint(ctypes.Structure):
    """This is a wrapper to C definition of a changepoint."""

    _fields_ = [
        ("significance_std", ctypes.c_double),
        ("changepoint", ctypes.c_size_t),
        ("triggertime", ctypes.c_size_t),
    ]


class _Changepoints(ctypes.Structure):
    """Hosts four changepoints, one per detectors"""

    _fields_ = [("changepoints", _Changepoint * NUMDETECTOR)]


class BftCWrapper(TriggerAlgorithm):
    """A wrapper to the C implementation of the BFT."""

    NDPOINTER = np.ctypeslib.ndpointer(dtype=ctypes.c_int64, ndim=2)

    def __str__(self):
        return ":fire:cBFT:fire:"

    def __init__(
        self,
        threshold_std: float,
        alpha: float,
        m: int,
        sleep: int,
        mu_min: float = 1.0,
        majority: int = 3,
    ):
        self.check_init_parameters(threshold_std, mu_min, alpha, m, sleep, majority)
        self.threshold_std = threshold_std
        self.mu_min = mu_min
        self.alpha = alpha
        self.m = m
        self.sleep = sleep
        self.majority = majority
        self._call = self.bind_bft_interface()

    def __call__(
        self,
        xss: npt.NDArray[np.int64],  # shape (4, _)
    ) -> Changepoint:
        cs = _Changepoints()
        _, xs_length = xss.shape
        error_code = self._call(
            ctypes.byref(cs),
            # TODO: avoid this fucking casting
            xss.astype(ctypes.c_int64),
            xs_length,
            self.threshold_std,
            self.mu_min,
            self.alpha,
            self.m,
            self.sleep,
            self.majority,
        )

        match error_code:
            case _Errors.INVALID_INPUT:
                raise ValueError("The inputs contain invalid entries.")
            case _Errors.INVALID_ALLOCATION:
                raise BufferError("Can't allocate memory for the algorithm.")

        significance_std = max([c.significance_std for c in cs.changepoints])
        changepoint = min([c.changepoint for c in cs.changepoints])
        triggertime = max([c.triggertime for c in cs.changepoints])
        return significance_std, changepoint, triggertime

    @staticmethod
    def check_init_parameters(
        threshold_std: float,
        mu_min: float,
        alpha: float,
        m: int,
        sleep: int,
        majority: int,
    ):
        """Checks validity of initialization arguments."""
        check = BftCWrapper.bind_bft_check_init_parameters()
        error_code = check(
            threshold_std,
            mu_min,
            alpha,
            m,
            sleep,
            majority,
        )
        match error_code:
            case _Errors.INVALID_INPUT:
                raise ValueError("The inputs contain invalid entries.")
        return

    @staticmethod
    def bind_bft_check_init_parameters() -> Callable:
        """Ctypes binding for C library interface."""
        bft_check_init_parameters = clib_bft.bft_check_init_parameters
        bft_check_init_parameters.restype = _Errors
        bft_check_init_parameters.argtypes = [
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_int,
        ]
        return bft_check_init_parameters

    def bind_bft_interface(self) -> Callable:
        """Ctypes binding for C library interface."""
        bft_interface = clib_bft.bft_interface
        bft_interface.restype = _Errors
        bft_interface.argtypes = [
            ctypes.POINTER(_Changepoints),
            self.NDPOINTER,
            ctypes.c_size_t,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_int,
        ]
        return bft_interface
