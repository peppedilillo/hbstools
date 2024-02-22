"""
An implementation of Poisson FOCuS with optimizations.
See Ward 2023 & Dilillo 2024.
"""

from math import inf
from math import log
from math import sqrt


class Curve:
    def __init__(self, x: float, b: float, t: int, m: float):
        self.x = x
        self.b = b
        self.t = t
        self.m = m


def ymax(curve, acc):
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


class PoissonFocus:
    """
    An implementation of Poisson-FOCuS with optimizations,
    see Ward et al., 2023. or Dilillo et al., 2024
    """

    def __init__(
        self,
        threshold_std: float,
        mu_min: float = 1.0,
    ):
        """
        Args:
            threshold_std: threshold value in standard deviations.
            mu_min: mumin value.
        """
        if mu_min < 1:
            raise ValueError("mumin must be greater or equal 1.0")
        if threshold_std <= 0:
            raise ValueError("threshold must be greater than 0.")

        self.ab_crit = 1 if mu_min == 1 else (mu_min - 1) / log(mu_min)
        self.threshold_llr = threshold_std**2 / 2
        self.global_max = 0.0
        self.time_offset = 0
        self.curve_list = []
        self.curve_list.append(Curve(inf, 0.0, 0, 0.0))
        self.curve_list.append(Curve(0, 0.0, 0, 0.0))

    def __call__(
        self,
        xs: list[int],
        bs: list[float],
    ):
        """
        Args:
            xs: a list of count data
            bs: a list of background values

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
