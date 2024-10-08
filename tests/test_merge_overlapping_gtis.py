import unittest

from hbstools.types import GTI

# fmt: off
EPSILON = 10**-7
inputs = [
    ([(1, 2)], EPSILON),
    ([(1, 2), (2, 3)], EPSILON),
    ([(1, 2), (3, 4)], EPSILON),
    ([(1, 2), (2, 3), (3, 4)], EPSILON),
    ([(1, 2), (2, 3), (4, 5)], EPSILON),
    ([(1, 2), (3, 4), (4, 5)], EPSILON),
    ([(1, 2), (3, 4), (5, 6)], EPSILON),
    ([(0, 5400), (5100, 10800), (10800, 16140), (16200, 21600), (21700, 27000)], EPSILON),
    ([(0, 5400), (5100, 7900), (8300, 10800), (10800, 13400), (13600, 16140),
      (16200, 18800), (19000, 21600), (21700, 24200), (24400, 27000),], EPSILON),
]

outputs = [
    [(1, 2)],
    [(1, 3)],
    [(1, 2), (3, 4)],
    [(1, 4)],
    [(1, 3), (4, 5)],
    [(1, 2), (3, 5)],
    [(1, 2), (3, 4), (5, 6)],
    [(0, 16140), (16200, 21600), (21700, 27000)],
    [(0, 7900), (8300, 13400), (13600, 16140), (16200, 18800), (19000, 21600), (21700, 24200), (24400, 27000)],
]
# fmt: on


def merge_overlapping_gtis(gtis: list[GTI], tolerance: float = 1.0) -> list[tuple]:
    """Merges CDI if they overlap. eg:
    F([(1, 2), (3, 4), (4, 5)]) = [(1, 2), (3, 5)]"""
    from math import isclose

    def overlap(x: GTI, y: GTI, abs_tol: float):
        """`x` overlaps `y` if `x` starts before or at the end of `y`"""
        assert x.start < y.stop
        return isclose(x.stop, y.start, abs_tol=abs_tol) or (y.start < x.stop)

    def f(xs: list[GTI]):
        if len(xs) == 2:
            first, second = xs
            if overlap(first, second, tolerance):
                return [GTI(first.start, second.stop)]
        if len(xs) > 2:
            first, second, *cdr = xs
            if overlap(first, second, tolerance):
                return f([GTI(first.start, second.stop)] + cdr)
            else:
                return [first] + f([second, *cdr])
        return xs

    return f(gtis)


class TestMergeOverlappingGTIs(unittest.TestCase):
    def test_equal(self):
        for x, y in zip(inputs, outputs):
            tuples, tolerance = x
            gtis = list(map(lambda g: GTI(*g), tuples))
            self.assertEqual(merge_overlapping_gtis(gtis, tolerance), y)


if __name__ == "__main__":
    unittest.main()
