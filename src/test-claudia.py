
# **** Example unittest tests
#      :LOGBOOK:
#      - Note taken on [2011-08-23 Tue 11:02] \\
#        First version is a straight port of the nose tests I already had
#      :END:
#      :PROPERTIES:
#      :tangle:   ../src/test-claudia.py
#      :END:

# #+srcname: unittest-claudia

import unittest
from claudia import CloudyModel

class ClaudiaTest(unittest.TestCase):
    def setUp(self):
        "set up test fixtures"
        CloudyModel.indir = '../testdata'
        self.model = CloudyModel('sample01')

    # def teardown_func():
    #     "tear down test fixtures"

    def test_doomed_to_fail():
        self.assertEquals(1, 2)

    def test_infilepath():
        self.assertEquals(self.model.infilepath, '../testdata/sample01.in')
