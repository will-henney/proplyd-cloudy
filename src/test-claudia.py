
# **** Example tests
#      :PROPERTIES:
#      :tangle:   ../src/test-claudia.py
#      :END:
# #+srcname: test-claudia-examples

import nose
from nose.tools import with_setup

def setup_func():
    "set up test fixtures"

def teardown_func():
    "tear down test fixtures"

@with_setup(setup_func, teardown_func)
def test():
    "test destined to fail"
    assert False
