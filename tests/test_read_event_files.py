import unittest

from hbstools.io import read_event_files


dataset_directory = "./data_100s_stronganomaly60s/"
df = read_event_files(dataset_directory)


class TestReadEvents(unittest.TestCase):
    def test_dataset_columns(self):
        self.assertTrue(len(df.columns) == 3)
        self.assertTrue("TIME" in df and "ENERGY" in df and "QUADID" in df)


if __name__ == "__main__":
    unittest.main()
