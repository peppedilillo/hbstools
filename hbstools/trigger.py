from typing import Callable, Type

import numpy as np
import pandas as pd

from hbstools.data import histogram
from hbstools.data import histogram_quadrants
import hbstools.triggers.bft as bft
import hbstools.triggers.bft_cwrap as bftc
import hbstools.triggers.poissonfocusdes as pfd
import hbstools.triggers.poissonfocusses_cwrap as pfsc
from hbstools.types import Changepoint
from hbstools.types import ChangepointMET
from hbstools.types import GTI


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
        case _:
            raise ValueError(f"Cannot find a binner for {algorithm}.")


def trigger_match(
    algorithm_params: dict,
) -> Type[
    pfd.PoissonFocusDes | bft.Bft | pfsc.PoissonFocusSesCwrapper | bftc.BftCWrapper
]:
    """Given a dictionary of parameters tries to find a suitable algoritm,
    giving precedence to algorithms with C implementations."""
    match algorithm_params:
        case {
            "threshold_std": _,
            "alpha": _,
            "beta": _,
            "m": _,
            "sleep": _,
            "mu_min": _,
            "majority": _,
            "t_max": _,
        }:
            return bft.Bft
        case {
            "threshold_std": _,
            "alpha": _,
            "beta": _,
            "m": _,
            "sleep": _,
            "mu_min": _,
            "t_max": _,
        }:
            return pfd.PoissonFocusDes
        case {
            "threshold_std": _,
            "alpha": _,
            "m": _,
            "sleep": _,
            "mu_min": _,
            "majority": _,
        }:
            return bftc.BftCWrapper
        case {
            "threshold_std": _,
            "alpha": _,
            "m": _,
            "sleep": _,
            "mu_min": _,
        }:
            return pfsc.PoissonFocusSesCwrapper
        case _:
            raise ValueError(
                "Cannot identify a fitting algorithm.\n" 
                "Check your configuration!"
            )


def trigger_counts_run(
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


def trigger_df_set(
    binning: float,
    skip: int,
    algorithm_params: dict,
) -> Callable:
    """Identifies a suitable algorithm and binner."""
    def trigger_df_run(data: pd.DataFrame, gti: GTI) -> list[ChangepointMET]:
        """Run the algorithm, restarting each time a trigger is found."""
        counts, bins = binner(data, gti, binning)
        # our choice of lambda implies the algorithm will be reset each time run_on_segment
        # calls the lambda function.
        cps = trigger_counts_run(lambda: algorithm(**algorithm_params), counts, bins, skip)
        # maps bin-steps to MET
        return list(map(lambda cp: (cp[0], bins[cp[1]], bins[cp[2]]), cps))

    algorithm = trigger_match(algorithm_params)
    # TODO: define a binner as an histogram with defined binning
    binner = _find_suitable_binner(algorithm)
    return trigger_df_run
