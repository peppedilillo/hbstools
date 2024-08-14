import unittest

from hbstools.format import format_results
from hbstools.types import GTI

binning = 1
alpha = 0.1
skip = 10
m = 5


format_params = {
    "intervals_duration_seconds": binning / alpha,
    "preinterval_ends_seconds": binning * m,
    "postinterval_start_seconds": binning * skip,
}


inputs = [
    {GTI(0, 50): [(0, 20.0, 25.0)]},
    {GTI(20, 80): [(0, 40.0, 42.0)]},
    {GTI(0, 50): [(0, 6.0, 7.0)]},
    {GTI(0, 50): [(0, 47.0, 48.0)]},
]

outputs = [
    [(10.0, 20.0, 20.0, 35.0, 35.0, 45.0)],
    #        ^^^ m = 5 seconds, changepoint = 20 s, trigger = 25 s
    [(27.0, 37.0, 40.0, 52.0, 52.0, 62.0)],
    [(0.0, 2.0, 6.0, 17.0, 17.0, 27.0)],  # left squashing
    [(33.0, 43.0, 47.0, 48.0, 48.0, 50.0)],  # right squashing
]


class TestFormatResult(unittest.TestCase):
    def test_equal(self):
        for x, y in zip(inputs, outputs):
            self.assertEqual(format_results(x, **format_params), y)


if __name__ == "__main__":
    unittest.main()
