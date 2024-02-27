import pandas as pd

import hbstools.triggers.poissonfocusdes as pfd
import hbstools.triggers.bft as bft
import hbstools.triggers.poissonfocusses_cwrap as pfsc
import hbstools.triggers.bft_cwrap as bftc
from hbstools.types import ChangepointMET, Changepoint, GTI
from hbstools.data import histogram, histogram_quadrants

from typing import Type, Callable
import numpy as np


def _find_suitable_binner(
    algorithm: Type[
        pfd.PoissonFocusDes | bft.Bft | pfsc.PoissonFocusSesCwrapper | bftc.BftCWrapper,
    ],
) -> Callable:
    """Given an algorithm type, finds a suitable binner"""
    match algorithm:
        case bftc.BftCWrapper | bft.Bft:
            return histogram_quadrants
        case pfsc.PoissonFocusSesCwrapper | pfd.PoissonFocusDes:
            return histogram


def get_algorithm(
    algorithm_params: dict,
) -> Type[
    pfd.PoissonFocusDes | bft.Bft | pfsc.PoissonFocusSesCwrapper | bftc.BftCWrapper
]:
    """Given a dictionary of parameters tries to find a suitable algoritm,
    giving precedence to algorithms with C implementations."""
    match algorithm_params:
        case {
            "threshold_std": _,
            "mu_min": _,
            "alpha": _,
            "beta": _,
            "m": _,
            "t_max": _,
            "sleep": _,
            "majority": _,
        }:
            return bft.Bft
        case {
            "threshold_std": _,
            "mu_min": _,
            "alpha": _,
            "beta": _,
            "m": _,
            "t_max": _,
            "sleep": _,
        }:
            return pfd.PoissonFocusDes
        case {
            "threshold_std": _,
            "mu_min": _,
            "alpha": _,
            "m": _,
            "sleep": _,
            "majority": _,
        }:
            return bftc.BftCWrapper


def _run_on_segment(
    init_algorithm: Callable,
    counts: np.ndarray,
    bins: np.ndarray,
    skip: int,
) -> list[Changepoint]:
    """Runs on binned data restarting the algorithm after skip bin-steps."""

    def helper(cs, bs, skip, acc):
        """Recursion helper"""
        match cs.shape:
            case (4, 0):
                return []
            case (4, _):
                # NOTE: we call init_algorithm() before launching it on the counts `cs`.
                # The rationale behind this design is to leave the user the freedom to
                # decide if some code (e.g., an initializer, a reset signal, whatever)
                # must run before the algorithm is actually launched on the data.
                s, cp, tt = init_algorithm()(cs)
                t = [(s, acc + cp, acc + tt)] if tt >= cp else []
                return t + helper(
                    cs[:, tt + skip :], bs[tt + skip :], skip, acc + tt + skip
                )
            case (0,):
                return []
            case (_,):
                s, cp, tt = init_algorithm()(cs)
                t = [(s, acc + cp, acc + tt)] if tt >= cp else []
                return t + helper(
                    cs[tt + skip :], bs[tt + skip :], skip, acc + tt + skip
                )
            case _:
                raise ValueError("Cannot match shape of counts with anything.")

    return helper(counts, bins, skip, 0)


def set(
    algorithm_params: dict,
) -> Callable:
    """Prepares the algorithm and runs on data.

    The outer layer identifies a suitable algorithm and binner.
    """

    def run(data: pd.DataFrame, gti: GTI, binning, skip: int) -> list[ChangepointMET]:
        """Run the algorithm, restarting each time a trigger is found."""
        counts, bins = binner(data, gti, binning)
        # our choice of lambda implies the algorithm will be reset each time run_on_segment
        # calls the lambda function.
        cps = _run_on_segment(lambda: algorithm(**algorithm_params), counts, bins, skip)
        # maps bin-steps to MET
        return list(map(lambda cp: (cp[0], bins[cp[1]], bins[cp[2]]), cps))

    algorithm = get_algorithm(algorithm_params)
    binner = _find_suitable_binner(algorithm)
    return run
