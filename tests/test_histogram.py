from math import isclose
import unittest

import numpy as np
import pandas as pd

from hbstools.data import histogram
from hbstools.types import GTI


inputs = [
    (pd.DataFrame({"TIME": [0.05,0.1, 0.2,0.35]}), GTI(0.05, 0.28), 0.1),
    (pd.DataFrame({"TIME": [0.06, 0.3]}), GTI(0.05, 0.28), 0.1),
]

outputs = [
    [np.array([2, 1, 1]), np.array([0.05, 0.15, 0.25, 0.35])],
    [np.array([1, 0, 1]), np.array([0.05, 0.15, 0.25, 0.35])],
]


class TestMakeBins(unittest.TestCase):
    def test_len_bins(self):
        for x, y in zip(inputs, outputs):
            expected_counts, expected_bins = y
            counts, bins = histogram(*x)

            self.assertEqual(len(expected_bins), len(bins))

    def test_len_counts(self):
        for x, y in zip(inputs, outputs):
            expected_counts, expected_bins = y
            counts, bins = histogram(*x)

            self.assertEqual(len(expected_counts), len(counts))

    def test_bin_values(self):
        for x, y in zip(inputs, outputs):
            expected_counts, expected_bins = y
            counts, bins = histogram(*x)

            for i in range(len(bins)):
                self.assertTrue(isclose(bins[i], expected_bins[i]))

    def test_counts_values(self):
        for x, y in zip(inputs, outputs):
            expected_counts, expected_bins = y
            counts, bins = histogram(*x)

            for i in range(len(counts)):
                self.assertTrue(counts[i] == expected_counts[i])

    def test_counts_shorter_than_bins(self):
        for x, y in zip(inputs, outputs):
            counts, bins = histogram(*x)

            self.assertTrue(len(counts) == len(bins) - 1)

    def test_bins_equispaced(self):
        for x, y in zip(inputs, outputs):
            counts, bins = histogram(*x)

            self.assertTrue(all(map(lambda x: np.isclose(x, 0), np.diff(bins) - 0.1)))


if __name__ == "__main__":
    unittest.main()
