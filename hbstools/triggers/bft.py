"""
Wraps a number of independent Poisson-FOCuS algorithms with DES background
estimate and triggers if more than half of them are over threshold.
"""

from typing import Callable

import numpy.typing as npt

from hbstools.triggers import TriggerAlgorithm
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

    DETECTOR_NUMBER = 4

    def __str__(self):
        return "BFT"

    def __init__(
        self,
        thr_std: float,
        alpha: float,
        beta: float,
        m: int,
        sleep: int,
        majority: int,
        mu_min: float = 1.0,
        t_max: int | None = None,
        s_0: float | None = None,
        b_0: float | None = None,
    ):
        """
        Args:
            thr_std:  In standard deviation units.
            Must be greater than 0.
            alpha: DES alpha (value) parameter.
            Must be greater than 0.
            beta: DES beta (slope) parameter.
            Must be non-negative.
            m: background estimate delay and forecast length.
            Must be a positive integer.
            sleep: dead time for automated s_0 initialization.
            Must be a non-negative integer.
            majority: Sets minimum number of dets. over threshold for trigger.
            Must be comprised between 1 and self.DETECTOR_NUMBER.
            mu_min: FOCuS mu_min parameter.
            Must not be smaller than 1.0
            t_max:  Maximum changepoint duration, for quality control.
            Must be a positive integer.
            s_0: DES init count rate parameter.
            Must be greater than zero.
            b_0: DES init slope parameter.
            Must be a non-negative number, set to 0.

        Optional arguments are set off by default.
        """
        self.majority = majority
        self.quadrantmask = [1 for i in range(self.DETECTOR_NUMBER)]
        self.fs = [
            PoissonFocusDes(
                thr_std=thr_std,
                alpha=alpha,
                beta=beta,
                m=m,
                sleep=sleep,
                mu_min=mu_min,
                t_max=t_max,
                s_0=s_0,
                b_0=b_0,
            )
            for _ in range(self.DETECTOR_NUMBER)
        ]

    @fold_changepoints
    def __call__(
        self,
        xss: npt.NDArray,  # shape (4, _)
    ) -> list[Changepoint]:
        assert len(xss) == self.DETECTOR_NUMBER
        changes = [(0.0, 0), (0.0, 0), (0.0, 0), (0.0, 0)]
        t_length = len(xss[0])
        for t in range(t_length):
            changes = self.step(xss[:, t])
            if sum(self.quadrantmask) < self.majority:
                raise ValueError("Not enough working quadrants for majority.")
            elif len([*filter(lambda c: c[0] > 0, changes)]) >= self.majority:
                break
        return [*map(lambda c: (c[0], t - c[1] + 1, t), changes)]

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
        return [self._step(det_id, xts) for det_id in range(self.DETECTOR_NUMBER)]
