import unittest

import hbstools.triggers.bft as bft
from hbstools import search
from hbstools.data import catalog
from hbstools.trigger import trigger_match
from hbstools.types import GTI

data_paths = ["./data_100s_stronganomaly60s/"]

TRIGTIME = 60
gti = GTI(0.0, 100.0)
configuration = {
    "binning": 0.1,
    "energy_lims": (20, 300),
    "skip": 10,
    "algorithm_params": {
        "threshold_std": 4.5,
        "mu_min": 1.1,
        "alpha": 0.005,
        "beta": 0.0,
        "m": 40,
        "t_max": 40,
        "sleep": 120,
        "majority": 1,
    },
}


class TestBFT(unittest.TestCase):
    def test_is_found(self):
        algorithm = trigger_match(configuration["algorithm_params"])
        self.assertTrue(algorithm is bft.Bft)

    def test_it_runs(self):
        dataset = catalog(data_paths)
        results = search(dataset, configuration)

        self.assertTrue(len(results) == 1)
        self.assertTrue(abs(results["event_start"].iloc[0] - TRIGTIME) < 5)

    if __name__ == "__main__":
        unittest.main()
