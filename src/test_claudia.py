
# **** Example unittest tests
#      :LOGBOOK:
#      - Note taken on [2011-08-23 Tue 11:10] \\
#        Note that we had to use test_claudia.py not test-claudia.py since the latter is not a valid module name.
#      - Note taken on [2011-08-23 Tue 11:02] \\
#        First version is a straight port of the nose tests I already had
#      :END:
#      :PROPERTIES:
#      :tangle:   ../src/test_claudia.py
#      :END:

# #+srcname: unittest-claudia

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
