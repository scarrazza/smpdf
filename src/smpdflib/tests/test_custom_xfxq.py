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


class TestXFXQ(unittest.TestCase):

    def setUp(self):
        try:
            init_multiprocessing()
        except RuntimeError:
            pass
        logging.basicConfig(level=logging.DEBUG)

    def test_lincomb_xfxq(self):
         obs = make_observable('data/applgrid/atlas-incljets-eta7.root', 'NLO')
         pdf = PDF('1000rep')
         lincomb = np.eye(1000,2)
         lpdf = LincombPDF(pdf, lincomb)
         produce_results([lpdf], [obs])



if __name__ == '__main__':
    unittest.main()