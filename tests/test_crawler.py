from pathlib import Path
import unittest

from mercury.mercury import crawler

# fmt: off
inputs = [
    ("crawler_test", ["file1", "file2"]),
    ("crawler_test", ["file1", "file2", "file3"], 1),
    ("crawler_test", ["file1", "file2", "file3"], 2),
    ("crawler_test", ["file1", "file2", "file3"], 3),
]

outputs = [
    {
        (Path("crawler_test/d4/file1"), Path("crawler_test/d4/file2")),
        (Path("crawler_test/d5/file1"), Path("crawler_test/d5/file2")),
    },
    {
        (Path("crawler_test/d4/file1"), Path("crawler_test/d4/file2"), Path("crawler_test/d4/file3")),
    },
    {
        (Path("crawler_test/d4/file1"), Path("crawler_test/d4/file2"), Path("crawler_test/d4/file3")),
        (Path("crawler_test/d2/d2_3/file1"), Path("crawler_test/d2/d2_3/file2"), Path("crawler_test/d2/d2_3/file3")),
    },
    {
        (Path("crawler_test/d4/file1"), Path("crawler_test/d4/file2"), Path("crawler_test/d4/file3")),
        (Path("crawler_test/d2/d2_3/file1"), Path("crawler_test/d2/d2_3/file2"), Path("crawler_test/d2/d2_3/file3")),
        (Path("crawler_test/d1/d1_1/d1_1_1/file1"), Path("crawler_test/d1/d1_1/d1_1_1/file2"), Path("crawler_test/d1/d1_1/d1_1_1/file3")),
    },
]


class TestMercuryCrawler(unittest.TestCase):
    def test_results(self):
        for x, y in zip(inputs, outputs):
            results = crawler(*x)
            self.assertEqual(results, y)


if __name__ == "__main__":
    unittest.main()
