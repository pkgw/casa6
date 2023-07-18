import os, shutil
import unittest
from casatestutils import sparse_check

class test_data_download(unittest.TestCase):
    @classmethod
    def tearDownClass(cls):
        for dataset in ["ngc5921.ms", "pm_ngc5921.ms","gaincaltest2.ms","gaincaltest2.ms.G0","gaincaltest2.ms.T0"]:
            if os.path.exists(dataset):
                shutil.rmtree(dataset)

    def test_single_download_similar_name(self):
        sparse_check.download_data(["ngc5921.ms"])
        self.assertTrue(os.path.exists("ngc5921.ms"))
        self.assertFalse(os.path.exists("pm_ngc5921.ms"))

    def test_single_download_data_ms_no_table(self):
        sparse_check.download_data(["gaincaltest2.ms"])
        self.assertTrue(os.path.exists("gaincaltest2.ms"))
        self.assertFalse(os.path.exists("gaincaltest2.ms.G0"))
        self.assertFalse(os.path.exists("gaincaltest2.ms.T0"))

if __name__ == '__main__':
    unittest.main()
