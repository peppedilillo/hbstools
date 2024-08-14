import unittest

from hbstools.data import catalog
from hbstools.search import search
from hbstools.trigger import trigger_match
import hbstools.triggers.poissonfocusses_cwrap as cpfs
from hbstools.types import GTI

data_paths = ["./data_100s_stronganomaly60s/"]

TRIGTIME = 60
gti = GTI(0.0, 100.0)
configuration = {
    "binning": 0.1,
    "energy_lims": (20, 300),
    "skip": 10,
    "algorithm_params": {
        "threshold_std": 5.0,
        "alpha": 0.005,
        "m": 40,
        "sleep": 120,
        "mu_min": 1.1,
    },
}


class TestCPFS(unittest.TestCase):
    def test_is_found(self):
        algorithm = trigger_match(configuration["algorithm_params"])
        self.assertTrue(algorithm is cpfs.PoissonFocusSesCwrapper)

    def test_it_runs(self):
        dataset = catalog(data_paths)
        results = search(dataset, configuration)

        self.assertTrue(len(results) == 1)
        self.assertTrue(abs(results[0].start - TRIGTIME) < 5)


if __name__ == "__main__":
    unittest.main()
