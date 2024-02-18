from math import isclose
import unittest

import numpy as np

from hbstools import Search

inputs = [
    (0.05, 0.28, 0.1),
]

outputs = [
    np.array([0.05, 0.15, 0.25, 0.35]),
]


class TestMakeBins(unittest.TestCase):
    def test_equal(self):
        for x, y in zip(inputs, outputs):
            result = Search.make_bins(*x)
            self.assertEqual(len(result), len(y))
            for i in range(len(result)):
                self.assertTrue(isclose(result[i], y[i]))


if __name__ == "__main__":
    unittest.main()
