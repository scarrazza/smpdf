# -*- coding: utf-8 -*-
"""
Created on Wed May 20 13:21:37 2015

@author: zah
"""

import unittest


from smpdflib.core import PDF, Observable, Result


class TestObs(unittest.TestCase):

    def test_prop(self):
         obs = Observable('data/applgrid/atlas-incljets-eta7.root', 1)
         mQ = [1,2,3,4,5,6]
         obs._meanQ = mQ
         self.assertEquals(list(obs.meanQ), mQ)


if __name__ == '__main__':
    unittest.main()