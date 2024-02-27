"""
Wraps a number of independent Poisson-FOCuS algorithms with DES background
estimate and triggers if more than half of them are over threshold.
"""

from typing import Callable

import numpy.typing as npt

from hbstools.triggers.poissonfocus import PoissonFocus
from hbstools.triggers.poissonfocusdes import PoissonFocusDes
from hbstools.types import Change
from hbstools.types import Changepoint

_QUADRANTS_NUMBER = 4


def fold_changepoints(func: Callable):
    """A decorator for reducing a list of changepoints from different detectors
    to a single changepoint"""

    def wrapper(*args, **kwargs) -> Changepoint:
        """Helper function"""
        changepoints = func(*args, **kwargs)
        significance_std = max([c[0] for c in changepoints])
        changepoint = min([c[1] for c in changepoints])
        triggertime = max([c[2] for c in changepoints])
        return significance_std, changepoint, triggertime

    return wrapper


class Bft:
    """BFT stands for Big Focus Trigger :^). It is a manager of multiple,
    independent FOCuS algorithms, with autonoumous background estimate by
    double exponential smoothing."""

    def __init__(
        self,
        threshold_std: float,
        alpha: float,
        beta: float,
        m: int,
        sleep: int,
        mu_min: float = 1.0,
        t_max: int | None = None,
        s_0: float | None = None,
        b_0: float | None = None,
        majority: int = 3,
    ):
        self.check_init_parameters(
            threshold_std, alpha, beta, m, sleep, t_max, mu_min, s_0, b_0, majority
        )
        self.majority = majority
        self.fs = [
            PoissonFocusDes(
                threshold_std=threshold_std,
                alpha=alpha,
                beta=beta,
                m=m,
                sleep=sleep,
                mu_min=mu_min,
                t_max=t_max,
                s_0=s_0,
                b_0=b_0,
            )
            for _ in range(_QUADRANTS_NUMBER)
        ]

    @fold_changepoints
    def __call__(
        self,
        xss: npt.NDArray,  # shape (4, _)
    ) -> list[Changepoint]:
        det_ids, _ = xss.shape
        changes = [(0.0, 0), (0.0, 0), (0.0, 0), (0.0, 0)]
        t_length = len(xss[0])
        for t in range(t_length):
            changes = self.step([xss[det_id, t] for det_id in range(det_ids)])
            if len([*filter(lambda c: c[0] > 0, changes)]) >= self.majority:
                break
        return [*map(lambda c: (c[0], t - c[1] + 1, t), changes)]

    @staticmethod
    def check_init_parameters(
        threshold_std, alpha, beta, m, sleep, t_max, mu_min, s_0, b_0, majority
    ):
        """Checks validity of initialization arguments."""
        PoissonFocus.check_init_parameters(
            threshold_std,
            mu_min,
        )
        PoissonFocusDes.check_init_parameters(
            threshold_std, alpha, beta, m, sleep, mu_min, t_max, s_0, b_0
        )
        if majority < 1 or majority > _QUADRANTS_NUMBER:
            raise ValueError(
                "majority should be comprised between 1 and `detector_number`"
            )
        return

    def step(self, xts: list[int]) -> list[Change]:
        """Basic algorithm step, i.e. call subalgorithms and asks them to do
        their thing."""
        # returns 4 quality-checked changes
        return [f.step(x_t) for f, x_t in zip(self.fs, xts)]
