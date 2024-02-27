import unittest

from hbstools import Search
from hbstools.trigger import get_algorithm
import hbstools.triggers.bft_cwrap as bftc
from hbstools.types import GTI

dataset_directory = "./data_100s_stronganomaly60s/"

TRIGTIME = 60
gti = GTI(0.0, 100.0)
binning = 0.1
energy_lims = (20, 300)
skip = 10
algorithm_params = {
    "threshold_std": 4.5,
    "mu_min": 1.1,
    "alpha": 0.005,
    "m": 40,
    "sleep": 120,
    "majority": 1,
}


class TestCBFT(unittest.TestCase):
    def test_is_found(self):
        algorithm = get_algorithm(algorithm_params)
        self.assertTrue(algorithm is bftc.BftCWrapper)

    def test_it_runs(self):
        search = Search(binning, skip, energy_lims, algorithm_params)
        results = search([dataset_directory])

        self.assertTrue(len(results) == 1)
        self.assertTrue(abs(results["event_start"].iloc[0] - TRIGTIME) < 5)


if __name__ == "__main__":
    unittest.main()
