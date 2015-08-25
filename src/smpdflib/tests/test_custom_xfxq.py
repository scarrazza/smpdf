# -*- coding: utf-8 -*-
"""
Created on Wed May 20 13:21:37 2015

@author: zah
"""

import unittest
import logging

import numpy as np

from smpdflib.initialization import init_multiprocessing
from smpdflib.core import PDF, LincombPDF, produce_results, make_observable


class BadError(Exception): pass

class Bad(PDF):
    CUSTOM_XFXQ = True
    def __len__(self):
        return 1

class Bad1(Bad):

    def xfxQ(self, replica, fl, x, Q):
        raise BadError()

class Bad2(Bad):
    def xfxQ(self, replica, fl, x, Q):
        return "patata"


class TestXFXQ(unittest.TestCase):

    def setUp(self):
        try:
            init_multiprocessing()
        except RuntimeError:
            pass
        logging.basicConfig(level=logging.DEBUG)

    def test_bad_values(self):
        obs = make_observable('data/applgrid/atlas-incljets-eta7.root', 'NLO')
        bad1 = Bad1('1000rep')
        bad2 = Bad2('1000rep')
        with self.assertRaises(BadError):
            produce_results([bad1], [obs])
        with self.assertRaises(TypeError):
            produce_results([bad2], [obs])

    def test_lincomb_xfxq(self):
         obs = make_observable('data/applgrid/atlas-incljets-eta7.root', 'NLO')
         pdf = PDF('1000rep')
         lincomb = np.eye(len(pdf) - 1, 2)
         lpdf = LincombPDF(pdf, lincomb)
         lin, orig = produce_results([lpdf, pdf], [obs])
         self.assertTrue( np.all(orig._all_vals[[1,2]] == lin._all_vals))


if __name__ == '__main__':
    unittest.main()
