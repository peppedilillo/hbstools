"""
A wrapper to Poisson-FOCuS implementing background estimate via double
exponential smoothing.
"""

from collections import deque
from math import sqrt
from typing import Deque, Sequence

from hbstools.triggers import TriggerAlgorithm
from hbstools.triggers.poissonfocus import PoissonFocus
from hbstools.types import Change
from hbstools.types import Changepoint


class PoissonFocusDes(TriggerAlgorithm):
    """
    A wrapper to Poisson-FOCuS implementing background estimate via double
    exponential smoothing.
    """

    def __str__(self):
        return "PF+DES"

    def __init__(
        self,
        thr_std: float,
        alpha: float,
        beta: float,
        m: int,
        sleep: int,
        mu_min: float = 1.0,
        t_max: int | None = None,
        s_0: float | None = None,
        b_0: float | None = None,
    ):
        """
        An implementation of Poisson FOCuS with automatic background estimate
        through double exponential smoothing.

        Args:
            thr_std:  In standard deviation units.
            Must be greater than 0.
            alpha: DES alpha (value) parameter.
            Must be greater than 0.
            beta: DES beta (slope) parameter.
            Must be non-negative.
            m: background estimate delay and forecast length.
            Must be a positive integer.
            t_max:  Maximum changepoint duration, for quality control.
            Must be a positive integer, disabled by default.
            sleep: dead time for automated s_0 initialization.
            Must be a non-negative integer.
            mu_min: FOCuS mu_min parameter.
            Must not be smaller than 1.0
            s_0: DES init count rate parameter.
            Must be greater than zero.
            b_0: DES init slope parameter.
            Must be a non-negative number.

        Optional arguments are implied off by default.
        """
        self.focus = PoissonFocus(thr_std, mu_min)
        self.buffer: Deque = deque([], maxlen=m)
        self.s_t = None
        self.b_t = None
        self.lambda_t = 0.0
        self.m = m
        self.t_max = t_max
        self.sleep = m + 1 if sleep is None else sleep
        self.alpha = alpha
        self.beta = beta
        self.t = m + sleep
        self.s_0 = s_0
        self.b_0 = b_0
        self.schedule = "collect"

    def __call__(
        self,
        xs: Sequence[int],
    ) -> Changepoint:
        """
        Args:
            xs: a list of count data

        Returns:
            A 3-tuple: significance value (std. devs), changepoint,  and
            stopping iteration (trigger time).

        Raises:
            ValueError: if zero background is passed to the update function.
        """
        t = 0
        for t, x_t in enumerate(xs):
            significance, offset = self.step(x_t)
            if significance:
                return significance, t - offset + 1, t
        return 0.0, t + 1, t

    def des_initialize(self):
        """Initialize background estimate."""
        if self.s_0 is None:
            assert len(self.buffer) == self.m
            counts_sum = sum(self.buffer)
            self.s_t = counts_sum / self.m
            self.lambda_t = self.s_t
        else:
            self.s_t = self.s_0

        if self.b_0 is None:
            self.b_t = 0.0
        else:
            self.b_t = self.b_0
        return

    def des_update(self, x):
        """Updates background estimate."""
        s_t_1 = self.s_t
        b_t_1 = self.b_t
        self.s_t = self.alpha * x + (1 - self.alpha) * (s_t_1 + b_t_1)
        self.b_t = self.beta * (self.s_t - s_t_1) + (1 - self.beta) * b_t_1
        return self.s_t + self.m * self.b_t

    def qc(self, significance, offset) -> Change:
        """
        quality control.
        runs only if focus is over threshold and t_max is given.
        it recomputes the curve stack maximum, skipping curves earlier than
        t_max. useful with delayed background estimate.
        """
        if significance and offset < self.t_max:
            return sqrt(2 * significance), offset
        return 0.0, 0

    def step(self, x) -> Change:
        """Base algorithm step."""
        if self.schedule == "test":
            x_t_m = self.buffer.popleft()
            self.lambda_t = self.des_update(x_t_m)
            self.buffer.append(x)
            if self.lambda_t <= 0.0:
                raise ValueError("background should be positive")
            self.focus.update(x, self.lambda_t)
            return self.qc(self.focus.global_max, self.focus.time_offset)

        elif self.schedule == "update":
            x_t_m = self.buffer.popleft()
            self.lambda_t = self.des_update(x_t_m)
            self.buffer.append(x)

            self.t = self.t - 1
            if self.t == 0:
                self.schedule = "test"
            return 0.0, 0

        elif self.schedule == "collect":
            self.buffer.append(x)
            self.t = self.t - 1
            if self.t == self.sleep:
                assert len(self.buffer) == self.m
                self.des_initialize()
                self.schedule = "update" if self.sleep else "test"
            return 0.0, 0
        raise ValueError("Unknown task.")
