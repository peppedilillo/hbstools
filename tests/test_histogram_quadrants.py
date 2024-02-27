import unittest

from hbstools.data import get_data
from hbstools.data import histogram_quadrants

dataset_directory = "./data_100s_stronganomaly60s/"
generator = get_data([dataset_directory])


class TestHistogramQuadrants(unittest.TestCase):
    def test_counts_have_right_shape(self):
        for data, gti in generator:
            counts, bins = histogram_quadrants(data, gti, 0.1)
            self.assertTrue(counts.shape == (4, 1001))

    def test_still_get_counts_if_quadrant_is_not_working(self):
        for data, gti in generator:
            data = data[data["QUADID"] != 0]
            counts, bins = histogram_quadrants(data, gti, 0.1)
            self.assertTrue(counts.shape == (4, 1001))


if __name__ == "__main__":
    unittest.main()
