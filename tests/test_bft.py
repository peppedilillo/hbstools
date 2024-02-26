import unittest

from hbstools.io import read_event_files
from hbstools.trigger import trigger_algorithm
from hbstools.types import GTI

dataset_directory = "./data_100s_stronganomaly60s/"

df = read_event_files(dataset_directory)
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
    "majority": 3,
}


class TestBFT(unittest.TestCase):
    def test_dataset_columns(self):
        result = trigger_algorithm(algorithm_params)()(df, gti, binning, skip)
        print(result)


if __name__ == "__main__":
    unittest.main()
