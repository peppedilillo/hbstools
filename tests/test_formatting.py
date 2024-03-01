import unittest

from hbstools.format import format_results
from hbstools.types import GTI

binning = 1
alpha = 0.1
skip = 10
m = 5


format_params = {
    "intervals_duration_seconds": binning / alpha,
    "preinterva_ends_seconds": binning / m,
    "postinterval_start_seconds": binning * skip,
}


inputs = [
    {GTI(0, 50): [(0, 20.0, 25.0)]},
    {GTI(20, 80): [(0, 40.0, 42.0)]},
    {GTI(0, 50): [(0, 6.0, 7.0)]},
    {GTI(0, 50): [(0, 47.0, 48.0)]},
]

outputs = [
    [(-10.0, 0.0, 20.0, 35.0, +15.0, +25.0)],
    #        ^^^ m = 5 seconds, changepoint = 20 s, trigger = 25 s
    [(-13.0, -3.0, 40.0, 52.0, +12.0, +22.0)],
    [(-6.0, -4.0, 6.0, 17.0, +11.0, +21.0)],
    [(-14.0, -4.0, 47.0, 48.0, +1.0, +3.0)],
]


class TestFormatResult(unittest.TestCase):
    def test_equal(self):
        for x, y in zip(inputs, outputs):
            self.assertEqual(format_results.__wrapped__(x, **format_params), y)


if __name__ == "__main__":
    unittest.main()
