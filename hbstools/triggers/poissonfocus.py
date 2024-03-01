"""
An implementation of Poisson FOCuS with optimizations.
See Ward 2023 & Dilillo 2024.
"""

from math import inf
from math import log
from math import sqrt
from typing import Sequence

from hbstools.triggers import TriggerAlgorithm


class Curve:
    """A Poisson-FOCuS curve supporting maximizer optimization."""

    def __init__(self, x: float, b: float, t: int, m: float):
        self.x = x
        self.b = b
        self.t = t
        self.m = m


def ymax(curve, acc):
    """Maximum of a curve, supporting accumulator optimization."""
    x = acc.x - curve.x
    b = acc.b - curve.b
    assert x > b
    return x * log(x / b) - (x - b)


def dominate(p, q, acc):
    """
    returns if p dominates q
    """
    area = (acc.x - p.x) * (acc.b - q.b) - (acc.x - q.x) * (acc.b - p.b)
    return +1 if area > 0 else -1


class PoissonFocus(TriggerAlgorithm):
    """
    An implementation of Poisson-FOCuS with optimizations,
    see Ward et al., 2023. or Dilillo et al., 2024
    """

    def __str__(self):
        return "PF"

    def __init__(
        self,
        threshold_std: float,
        mu_min: float = 1.0,
    ):
        """
        Args:
            threshold_std: Threshold value in standard deviations.
            Must be a positive number.
            mu_min: Kills old changepoints and keep amortized memory costant,
            at cost of a loss of sensitivity.
            Must be greater than 1.

        Optional arguments are implied off by default.
        """
        self.ab_crit = 1 if mu_min == 1.0 else (mu_min - 1.0) / log(mu_min)
        self.threshold_llr = threshold_std**2 / 2
        self.global_max = 0.0
        self.time_offset = 0
        self.curve_list = []
        self.curve_list.append(Curve(inf, 0.0, 0, 0.0))
        self.curve_list.append(Curve(0, 0.0, 0, 0.0))

    def __call__(
        self,
        xs: Sequence[int],
        bs: Sequence[float],
    ):
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
        for t, (x_t, b_t) in enumerate(zip(xs, bs)):
            if b_t <= 0:
                raise ValueError("background rate must be greater than zero.")
            self.update(x_t, b_t)
            if self.global_max > self.threshold_llr:
                return sqrt(2 * self.global_max), -self.time_offset + t + 1, t
        return 0.0, t + 1, t

    def update(self, x, b):
        """Poisson-FOCuS update step."""
        self.global_max = 0.0
        self.time_offset = 0

        p = self.curve_list.pop(-1)
        acc = Curve(p.x + x, p.b + b, p.t + 1, p.m)
        while dominate(p, self.curve_list[-1], acc) <= 0:
            p = self.curve_list.pop(-1)

        if (acc.x - p.x) > self.ab_crit * (acc.b - p.b):
            acc.m = p.m + ymax(p, acc)
            self.maximize(p, acc)
            self.curve_list.append(p)
            self.curve_list.append(acc)
        else:
            self.curve_list = self.curve_list[:1]
            self.curve_list.append(Curve(0, 0.0, 0, 0.0))
        return

    def maximize(self, p, acc):
        """Poisson-FOCuS maximization step."""
        m = acc.m - p.m
        i = len(self.curve_list)
        while m + p.m >= self.threshold_llr:
            if m >= self.threshold_llr:
                self.global_max = m
                self.time_offset = acc.t - p.t
                break
            i -= 1
            p = self.curve_list[i]
            m = ymax(p, acc)
        return
