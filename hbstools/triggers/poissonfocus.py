"""
A wrapper to Poisson-FOCuS implementing background estimate via double
exponential smoothing.
"""

from collections import deque
from math import inf
from math import log
from math import sqrt
from typing import Deque, Sequence


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
        for t, (x_t, b_t) in enumerate(zip(xs, bs)):
            if b_t <= 0:
                raise ValueError("background rate must be greater than zero.")
            self.update(x_t, b_t)
            if self.global_max > self.threshold_llr:
                return sqrt(2 * self.global_max), -self.time_offset + t + 1, t
        return 0.0, len(xs) + 1, len(xs)

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


class PoissonFocusDES:
    def __init__(
        self,
        threshold: float,
        alpha: float,
        beta: float,
        m: int,
        mu_min: float = 1.0,
        t_max: int | None = None,
        sleep: int | None = None,
        s_0: float | None = None,
        b_0: float | None = None,
    ):
        """
        An implementation of Poisson FOCuS with automatic background estimate
        through double exponential smoothing.

        Args:
            threshold:  in standard deviation units.
            alpha: DES alpha (value) parameter. Must be greater than 0.
            beta: DES beta (slope) parameter. Must be non-negative.
            m: background estimate delay and forecast length.
            t_max:  maximum changepoint duration, quality control.
            disabled by default.
            mu_min: FOCuS mu_min parameter. defaults to 1.
            sleep: dead time for automated s_0 initialization.
            if provided, it must be greater than m; else, defaults to m + 1.
            s_0: DES init count parameter.
            must be greater than 0.
            defaults to averaged over first `sleep - m` counts.
            b_0: DES init slope parameter. must be greater or equal than 0.
            defaults to 0.
        """
        if alpha < 0.0:
            raise ValueError("alpha must be non negative.")
        if beta < 0.0:
            raise ValueError("beta must be non negative.")
        if m < 0:
            raise ValueError("m must be non negative.")
        if (t_max is not None) and (t_max <= 0):
            raise ValueError("t_max must be greater than zero.")
        if (sleep is not None) and (sleep <= m):
            raise ValueError("sleep must be greater than m.")
        if (s_0 is not None) and (s_0 < 0):
            raise ValueError("s_0 must be non negative.")
        if (b_0 is not None) and (b_0 < 0):
            raise ValueError("b_0 must be non negative.")

        self.focus = PoissonFocus(threshold, mu_min=mu_min)
        self.buffer: Deque = deque([])
        self.s_t = None
        self.b_t = None
        self.lambda_t = None
        self.t = None
        self.m = m
        self.t_max = t_max
        self.sleep = m + 1 if sleep is None else sleep
        self.alpha = alpha
        self.beta = beta
        self.s_0 = s_0
        self.b_0 = b_0

    def des_initialize(self):
        if self.s_0 is None:
            counts_sum = sum([self.buffer[i] for i in range(self.sleep - self.m)])
            delta_t = self.sleep - self.m
            self.s_t = counts_sum / delta_t
        else:
            self.s_t = self.s_0
        if self.b_0 is None:
            self.b_t = 0.0
        else:
            self.b_t = self.b_0
        return

    def des_update(self, x):
        s_t_1 = self.s_t
        b_t_1 = self.b_t
        self.s_t = self.alpha * x + (1 - self.alpha) * (s_t_1 + b_t_1)
        self.b_t = self.beta * (self.s_t - s_t_1) + (1 - self.beta) * b_t_1
        return

    def qc(self):
        """
        quality control.
        runs only if focus is over threshold and t_max is given.
        it recomputes the curve stack maximum, skipping curves earlier than
        t_max. useful with delayed background estimate.
        """
        if self.t_max is None or self.focus.time_offset < self.t_max:
            return sqrt(2 * self.focus.global_max), self.focus.time_offset
        return 0.0, 0

    def step(self, x):
        self.t = 0 if self.t is None else self.t + 1
        self.buffer.append(x)
        if self.t <= self.sleep:
            if self.t != self.sleep:
                return 0.0, 0
            # initialize des
            if self.t == self.sleep:
                self.des_initialize()
                for i in range(self.sleep - self.m):
                    self.buffer.popleft()

        x_t_m = self.buffer.popleft()
        self.des_update(x_t_m)
        self.lambda_t = self.s_t + self.m * self.b_t
        if self.lambda_t <= 0:
            raise ValueError("background rate must be greater than zero.")
        self.focus.update(x, self.lambda_t)
        if self.focus.global_max:
            significance, offset = self.qc()
            return significance, offset
        return 0.0, 0


def init(**kwargs):
    """
    Initializes focus with double exponential smoothing and returns a callable
    trigger function. Arguments are passed to PoissonFocusDES.
    """

    # noinspection PyUnusedLocal
    def run(xs: Sequence[int], bins: Sequence[float] | None = None):
        """
        Args:
            xs: a sequence of count data
            bins: this is provided for compatibility with debugger utility.

        Returns:
            A 3-tuple: significance value (std. devs), changepoint,  and
            stopping iteration (trigger time).

        Raises:
            ValueError: if zero background is passed to the update function.
        """
        focus_des = PoissonFocusDES(**init_parameters)
        for t, x_t in enumerate(xs):
            significance, offset = focus_des.step(x_t)
            if significance:
                return significance, t - offset + 1, t
        return 0.0, len(xs) + 1, len(xs)

    init_parameters = kwargs
    return run


class Debugger:
    def __init__(self, **kwargs):
        self.init_parameters = kwargs
        self.log = {
            "timestamps": [],
            "ts": [],
            "xs": [],
            "bs": [],
        }

    def __call__(self, xs: Sequence[int], bins: Sequence[float]):
        """
        Args:
            xs: an iterable of count data
            bins: an iterable of the left edges of the bins

        Returns:
            A 3-tuple: significance value (std. devs), changepoint,  and
            stopping iteration (trigger time).

        Raises:
            ValueError: if zero background is passed to the update function.
        """
        focus_des = PoissonFocusDES(**self.init_parameters)
        for t, x_t in enumerate(xs):
            significance, offset = focus_des.step(x_t)
            self.log_step(bins[t], t, x_t, focus_des.lambda_t)
            if significance:
                return significance, t - offset + 1, t
        return 0.0, len(xs) + 1, len(xs)

    def log_step(self, bin_, t, x_t, b_t):
        self.log["timestamps"].append(bin_)
        self.log["ts"].append(t)
        self.log["xs"].append(x_t)
        self.log["bs"].append(b_t)

    def get_log(self):
        return self.log
