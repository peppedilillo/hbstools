from pathlib import Path
import unittest

from hbstools.search import Search
from hbstools.trigger import match_algorithm
import hbstools.triggers.poissonfocusses_cwrap as cpfs
from hbstools.types import GTI

dataset_directory = Path(__file__).parent / "data_100s_stronganomaly60s/"

TRIGTIME = 60
gti = GTI(0.0, 100.0)
binning = 0.1
energy_lims = (20, 300)
skip = 100
algorithm_params = {
    "threshold_std": 5.0,
    "alpha": 0.005,
    "m": 40,
    "sleep": 120,
    "mu_min": 1.1,
}


class TestCPFS(unittest.TestCase):
    def test_is_found(self):
        algorithm = match_algorithm(algorithm_params)
        self.assertTrue(algorithm is cpfs.PoissonFocusSesCwrapper)

    def test_it_runs(self):
        search = Search(binning, skip, energy_lims, algorithm_params)
        results = search([dataset_directory])

        self.assertTrue(search.algorithm_type is cpfs.PoissonFocusSesCwrapper)
        self.assertTrue(len(results) == 1)
        self.assertTrue(abs(results["event_start"].iloc[0] - TRIGTIME) < 5)


if __name__ == "__main__":
    unittest.main()
