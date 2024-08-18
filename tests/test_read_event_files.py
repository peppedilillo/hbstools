import unittest

from hbstools.io import read_event_files

dataset_directory = "./data_100s_stronganomaly60s/out_lv1_cl.evt"
df = read_event_files(dataset_directory)


class TestReadEvents(unittest.TestCase):
    def test_dataset_columns(self):
        self.assertTrue(len(df.columns) == 4)
        self.assertTrue(
            "TIME" in df and "ENERGY" in df and "QUADID" in df and "EVTYPE" in df
        )


if __name__ == "__main__":
    unittest.main()
