import unittest

from hbstools.data import catalog
from hbstools.data import stream
from hbstools.types import GTI

inputs = {
    "data_300s_constant/data_300s_constant1",
    "data_300s_constant/data_300s_constant2",
    "data_300s_constant/data_300s_constant3",
    "data_300s_constant/data_300s_constant4",
    "data_300s_constant/data_300s_constant5",
    "data_300s_constant/data_300s_constant6",
    "data_300s_constant/data_300s_constant7",
}

# key = merging tolerance
# value = merged gtis
out_gtis = {
    1.0: [
        GTI(0, 79),
        GTI(83, 133),
        GTI(137, 187),
        GTI(191, 241),
        GTI(245, 300),
    ],
    0.5: [
        GTI(0, 79),
        GTI(83, 133),
        GTI(137, 161),
        GTI(162, 187),
        GTI(191, 216),
        GTI(217, 241),
        GTI(245, 282),
        GTI(283, 300),
    ],
}


class TestDataset(unittest.TestCase):
    def test_length(self):
        self.assertTrue(len(catalog(inputs)) == 13)

    def test_file(self):
        filenames = set([filename for _, filename in catalog(inputs)])
        self.assertTrue(all([fn in inputs for fn in filenames]))

    def test_gtis(self):
        dataset = catalog(inputs)
        for abs_tol, expected_merged_gti in out_gtis.items():
            expected_merged_gti = out_gtis[abs_tol]
            datastream = stream(dataset, abs_tol=abs_tol)
            for i, (df, gti) in enumerate(datastream):
                self.assertTrue(expected_merged_gti[i] == gti)
                self.assertTrue(df["TIME"].max() < gti.end)
                self.assertTrue(df["TIME"].min() >= gti.start)
            self.assertTrue(len(expected_merged_gti) == i + 1)


if __name__ == "__main__":
    unittest.main()
