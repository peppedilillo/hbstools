import unittest

from hbstools import search
from hbstools.data import catalog
from hbstools.trigger import trigger_match
import hbstools.triggers.bft_cwrap as cbft
from hbstools.types import GTI

data_paths = [
    (
        "data_100s_stronganomaly60s/out_lv1_cl.evt",
        "data_100s_stronganomaly60s/gti.fits",
    )
]


TRIGTIME = 60
gti = GTI(0.0, 100.0)
configuration = {
    "binning": 0.1,
    "en_lims": (20, 300),
    "skip": 10,
    "algorithm_params": {
        "thr_std": 4.5,
        "mu_min": 1.1,
        "alpha": 0.005,
        "m": 40,
        "sleep": 120,
        "majority": 1,
    },
}


class TestCBFT(unittest.TestCase):
    def test_is_found(self):
        algorithm = trigger_match(configuration["algorithm_params"])
        self.assertTrue(algorithm is cbft.BftCWrapper)

    def test_it_runs(self):
        dataset = catalog(data_paths)
        results = search(dataset, configuration)

        self.assertTrue(len(results) == 1)
        self.assertTrue(abs(results[0].start - TRIGTIME) < 5)


if __name__ == "__main__":
    unittest.main()
