"""
Wraps a number of independent Poisson-FOCuS algorithms with DES background
estimate and triggers if more than half of them are over threshold.
"""

from typing import Callable

import numpy.typing as npt

from hbstools.triggers import TriggerAlgorithm
from hbstools.triggers.poissonfocus import PoissonFocus
from hbstools.triggers.poissonfocusdes import PoissonFocusDes
from hbstools.types import Change
from hbstools.types import Changepoint


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


class Bft(TriggerAlgorithm):
    """BFT stands for Big Focus Trigger :^). It is a manager of multiple,
    independent FOCuS algorithms, with autonoumous background estimate by
    double exponential smoothing."""

    QUADRANTS_NUMBER = 4

    def __str__(self):
        return "BFT"

    def __init__(
        self,
        threshold_std: float,
        alpha: float,
        beta: float,
        m: int,
        sleep: int,
        mu_min: float = 1.0,
        majority: int = 3,
        t_max: int | None = None,
        s_0: float | None = None,
        b_0: float | None = None,
    ):
        self.check_init_parameters(
            threshold_std, alpha, beta, m, sleep, t_max, mu_min, s_0, b_0, majority
        )
        self.majority = majority
        self.quadrantmask = [1 for i in range(self.QUADRANTS_NUMBER)]
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
            for _ in range(self.QUADRANTS_NUMBER)
        ]

    @fold_changepoints
    def __call__(
        self,
        xss: npt.NDArray,  # shape (4, _)
    ) -> list[Changepoint]:
        assert len(xss) == self.QUADRANTS_NUMBER
        changes = [(0.0, 0), (0.0, 0), (0.0, 0), (0.0, 0)]
        t_length = len(xss[0])
        for t in range(t_length):
            changes = self.step(xss[:, t])
            if sum(self.quadrantmask) < self.majority:
                raise ValueError("Not enough working quadrants for majority.")
            elif len([*filter(lambda c: c[0] > 0, changes)]) >= self.majority:
                break
        return [*map(lambda c: (c[0], t - c[1] + 1, t), changes)]

    def check_init_parameters(
        self,
        threshold_std,
        alpha,
        beta,
        m,
        sleep,
        mu_min,
        majority,
        t_max,
        s_0,
        b_0,
    ):
        """Checks validity of initialization arguments."""
        PoissonFocus.check_init_parameters(
            threshold_std,
            mu_min,
        )
        PoissonFocusDes.check_init_parameters(
            threshold_std, alpha, beta, m, sleep, mu_min, t_max, s_0, b_0
        )
        if majority < 1 or majority > self.QUADRANTS_NUMBER:
            raise ValueError("majority should be comprised between 1 and `detector_number`")
        return

    def _step(self, det_id: int, xts) -> Change:
        """Runs a single algorithms if it has been working so far,
         then annotates if an error occurs. Returns a trivial change"""
        if not self.quadrantmask[det_id]:
            return 0.0, 0
        try:
            c = self.fs[det_id].step(xts[det_id])
        except ValueError:
            self.quadrantmask[det_id] = 0
            return 0.0, 0
        return c

    def step(self, xts: list[int] | npt.NDArray) -> list[Change]:
        """Basic algorithm step, i.e. call subalgorithms and asks them to do
        their thing."""
        # returns 4 quality-checked changes
        return [self._step(det_id, xts) for det_id in range(self.QUADRANTS_NUMBER)]
