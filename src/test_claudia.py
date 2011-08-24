
import unittest
from claudia import CloudyModel

class ClaudiaTestSample01(unittest.TestCase):
    def setUp(self):
        "set up test fixtures"
        self.model = CloudyModel('sample01', 
                                 indir='../testdata', 
                                 outdir='../testdata',
                                 skipsaves=[])

    # def teardown_func():
    #     "tear down test fixtures"

    def test_doomed_to_fail(self):
        self.assertEquals(1, 2)

    def test_infilepath(self):
        self.assertEquals(self.model.infilepath, '../testdata/sample01.in')
