import unittest
from pathlib import Path

from mercury.mercury import crawler

inputs = [
    ("crawler_test", ["file1", "file2"]),
    ("crawler_test", ["file1", "file2", "file3"], 1),
    ("crawler_test", ["file1", "file2", "file3"], 2),
    ("crawler_test", ["file1", "file2", "file3"], 3),
]

outputs = [
    [Path("crawler_test").joinpath("d4"), Path("crawler_test").joinpath("d5")],
    [Path("crawler_test").joinpath("d4")],
    [Path("crawler_test").joinpath("d4"), Path("crawler_test").joinpath("d2/d2_3")],
    [
        Path("crawler_test").joinpath("d4"),
        Path("crawler_test").joinpath("d2/d2_3"),
        Path("crawler_test").joinpath("d1/d1_1/d1_1_1"),
    ],
]


class TestMercuryCrawler(unittest.TestCase):
    def test_results(self):
        for x, y in zip(inputs, outputs):
            results = crawler(*x)
            self.assertEqual(len(results), len(y))
            for r in results:
                self.assertIn(r, y)


if __name__ == "__main__":
    unittest.main()
