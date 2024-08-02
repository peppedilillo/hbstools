import unittest

from hbstools.data import catalog
from hbstools.types import GTI

input_paths = [
    {
        "data_300s_constant/test_dataset_light1",
        "data_300s_constant/test_dataset_light2",
        "data_300s_constant/test_dataset_light3",
        "data_300s_constant/test_dataset_light4",
        "data_300s_constant/test_dataset_light5",
        "data_300s_constant/test_dataset_light6",
        "data_300s_constant/test_dataset_light7",
    },
    {
        "data_100s_stronganomaly60s",
    },
]

output_gtis = [
    [
        [
            GTI(0, 54),
        ],
        [GTI(51, 79), GTI(83, 108)],
        [GTI(108.5, 133), GTI(137, 161)],
        [GTI(162, 187), GTI(191, 216)],
        [GTI(217, 241), GTI(245, 270)],
        [GTI(265, 280)],
        [GTI(280, 282), GTI(283, 285), GTI(285.01, 300)],
    ],
    [
        [GTI(0, 100)],
    ],
]


class TestDataset(unittest.TestCase):
    def test_file(self):
        for paths in input_paths:
            filenames = set([filename for _, filename in catalog(paths)])
            self.assertTrue(all([fn in paths for fn in filenames]))

    def test_gtis(self):
        for paths, gtis in zip(input_paths, output_gtis):
            flattened_gtis = [gti for file_gtis in gtis for gti in file_gtis]
            output = catalog(paths)
            expected = [gti for gti, _ in output]
            self.assertTrue(all([r == o for r, o in zip(flattened_gtis, expected)]))
            self.assertTrue(len(expected) == len(flattened_gtis))


if __name__ == "__main__":
    unittest.main()
