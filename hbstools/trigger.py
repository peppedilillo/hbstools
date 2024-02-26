import pandas as pd

import hbstools.triggers.poissonfocusdes as pfd
import hbstools.triggers.bft as bft
import hbstools.triggers.poissonfocusses_cwrap as pfsc
import hbstools.triggers.bft_cwrap as bftc
from hbstools.types import ChangepointMET, Changepoint, GTI
from hbstools.data import histogram, histogram_quadrants

from typing import Type, Callable
import numpy as np


def _find_suitable_algorithm(
        algorithm_params: dict
) -> Type[pfd.PoissonFocusDes | bft.Bft | pfsc.PoissonFocusSesCwrapper | bftc.BftCWrapper]:
    """Given a dictionary of parameters tries to find a suitable algoritm,
    giving precedence to algorithms with C implementations."""
    match algorithm_params:
        case {'threshold_std': _, 'mu_min': _, 'alpha': _, 'm': m, 'sleep': _, 'majority': _}:
            return bftc.BftCWrapper


def _find_suitable_binner(
        algorithm: Type[pfd.PoissonFocusDes | bft.Bft | pfsc.PoissonFocusSesCwrapper | bftc.BftCWrapper,],
) -> Callable:
    """Given an algorithm type, finds a suitable binner"""
    match algorithm:
        case bftc.BftCWrapper | bft.Bft:
            return histogram_quadrants
        case pfsc.PoissonFocusSesCwrapper | pfd.PoissonFocusDes:
            return histogram


def trigger_algorithm(algorithm_params: dict):
    """Prepares the algorithm and runs on data.
    Curried so that you can decide whether to keep using an already initialized
    algorithm or to initialize it again, without having to identify the right
    implementation and binner anew.

    The outer layer identifies a suitable algorithm and binner.
    """
    def run_on_segment(
        initialized_algorithm,
        counts: np.ndarray,
        bins: np.ndarray,
        skip: int,
    ) -> list[Changepoint]:
        """Runs on binned data restarting the algorithm after skip bin-steps."""

        def helper(cs, bs, skip, acc):
            match cs.shape:
                case (4, 0):
                    return []
                case (4, _):
                    s, cp, tt = initialized_algorithm(cs, bs)
                    t = [(s, acc + cp, acc + tt)] if tt >= cp else []
                    return t + helper(cs[:, tt + skip:], bs[tt + skip:], skip, acc + tt + skip)
                case (1, 0):
                    return []
                case (1, _):
                    s, cp, tt = initialized_algorithm(cs, bs)
                    t = [(s, acc + cp, acc + tt)] if tt >= cp else []
                    return t + helper(cs[tt + skip:], bs[tt + skip:], skip, acc + tt + skip)

        return helper(counts, bins, skip, 0)

    def initialize_algorithm():
        """Initialize the algorithm."""
        initialized_algorithm = algorithm(**algorithm_params)

        def run_algorithm(data: pd.DataFrame, gti: GTI, binning, skip: int) -> list[ChangepointMET]:
            """Run the algorithm."""
            counts, bins = binner(data, gti, binning)
            cps = run_on_segment(initialized_algorithm, counts, bins, skip)
            # maps bin-steps to MET
            return list(map(lambda cp: (cp[0], bins[cp[1]], bins[cp[2]]), cps))

        return run_algorithm

    algorithm = _find_suitable_algorithm(algorithm_params)
    binner = _find_suitable_binner(algorithm)
    return initialize_algorithm
