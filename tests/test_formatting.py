import unittest

from hbstools.search import Search
from hbstools.types import GTI

parameters = {
    "binning": 1,
    "alpha": 0.1,
    "m": 5,
    "skip": 10,
}

configuration = {
    "binning": 1,
    "skip": 10,
    "energy_lims": (None, None),
    "trigger_params": {
        "alpha": 0.1,
        "m": 5,
    },
}


inputs = [
    ((0, 20.0, 25.0), GTI(0, 50)),
    ((0, 40.0, 42.0), GTI(20, 80)),
    ((0, 6.0, 7.0), GTI(0, 50)),
    ((0, 47.0, 48.0), GTI(0, 50)),
]

outputs = [
    (-10.0, 0.0, 20.0, 35.0, +15.0, +25.0),
    (-13.0, -3.0, 40.0, 52.0, +12.0, +22.0),
    (-6.0, -4.0, 6.0, 17.0, +11.0, +21.0),
    (-14.0, -4.0, 47.0, 48.0, +1.0, +3.0),
]


class TestFormatResult(unittest.TestCase):
    def test_equal(self):
        search = Search(**configuration)
        for x, y in zip(inputs, outputs):
            print(search.format_result(*x))
            self.assertEqual(search.format_result(*x), y)


if __name__ == "__main__":
    unittest.main()
