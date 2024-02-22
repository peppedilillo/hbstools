"""
Wraps a number of independent Poisson-FOCuS algorithms with DES background
estimate and triggers if more then half of them are over threshold.
"""
from hbstools.triggers.poissonfocusdes import PoissonFocusDES


class Bft:
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
        detectors_number: int = 4,
    ):
        if detectors_number < 0:
            raise ValueError("detectors number must be a positive integer.")

        self.detectors_number = detectors_number
        self.focusexp_params = {
            "threshold_std": threshold_std,
            "alpha": alpha,
            "beta": beta,
            "m": m,
            "sleep": sleep,
            "mu_min": mu_min,
            "t_max": t_max,
            "s_0": s_0,
            "b_0": b_0,
        }
        self.fs = [
            PoissonFocusDES(**self.focusexp_params)
            for _ in range(self.detectors_number)
        ]

    def __call__(
        self,
        xss: Sequence[Sequence[int]],
        bins: Sequence[float] | None = None,
    ) -> list[tuple]:
        if len(xss) != self.detectors_number:
            raise ValueError(
                "the number of time series to analyze must match detector number"
            )
        if not reduce(lambda acc, x: acc and len(x) == len(xss[0]), xss, True):
            raise ValueError("time series must have same length")

        t_length = len(xss[0])
        for t in range(t_length):
            changes = self.step([xss[i][t] for i in range(self.detectors_number)])

            if len(changes) > self.detectors_number // 2:
                break
        return [
            (lambda s, o, t_: (s, t_ - o + 1, t_))(*changes[det_num], t)
            if (det_num in changes)
            else (0, t + 1, t)
            for det_num in range(self.detectors_number)
        ]

    @staticmethod
    def qc(changes: list[tuple]) -> dict[tuple]:
        return {
            det_num: (significance, offset)
            for det_num, (significance, offset) in enumerate(changes)
            if significance > 0
        }

    def step(self, xts: list[int]) -> dict[tuple]:
        if len(xts) != self.detectors_number:
            raise ValueError(
                "the number of data points to analyze must be equal to the number of detectors"
            )

        # returns 4 quality-checked changes
        return self.qc([f.step(x_t) for f, x_t in zip(self.fs, xts)])
