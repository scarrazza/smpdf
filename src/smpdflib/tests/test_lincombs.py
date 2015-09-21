# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 18:01:18 2015

@author: zah
"""

from collections import OrderedDict
import unittest

import numpy as np

from smpdflib.reducedset import merge_lincombs

class TestLincombs(unittest.TestCase):
    def test_merge_lincombs(self):
        desc1 = OrderedDict(( 
                            ('a', OrderedDict(( 
                                              (1,2),
                                              (2,4),
                                              (3,6),
                                              (4,6),
                                              (5,7))
                                  )),
                            ('b', OrderedDict(( 
                                              (1,7),
                                              (2,7),
                                              (3,10),
                                              (4,12),
                                              (5,14))
                          )),
                        
                            ))

        desc2 = OrderedDict(( 
                            ('a', OrderedDict(( 
                                              (2,2),
                                              (3,5),
                                              (4,5),
                                  ))),
                            ('b', OrderedDict(( 
                                              (1,6),
                                              (2,7),
                          ))),
                        
                            ))
        nrep = 2
        l1 = np.ones((nrep, 14))
        l2 = np.zeros((nrep, 7))
        final_lincomb, final_desc = merge_lincombs(l1,l2, desc1,desc2)
        self.assertEqual(final_desc, OrderedDict([
                        ('a', OrderedDict([(1, 2), 
                                           (2, 6), (3, 11), 
                                           (4, 11), (5, 12)])), 
                        ('b', OrderedDict([(1, 13), (2, 14), (3, 17), 
                                           (4, 19), (5, 21)]))])
                        )
        self.assertTrue((final_lincomb == np.array([[ 1,  1,  1,  1,  0,  0, 1,
                                                    1,  0,  0,  0,  1,  0,  0, 
                                                    1,  1,  1,  1,
                                                    1,  1,  1,],
                                                   [ 1,  1,  1,  1,  0,  0,  1,
                                                    1,  0,  0,  0,  1,  0,  0, 
                                                    1,  1,  1,  1, 1,  1,  1,]]
                                                    )).all())

if __name__ == '__main__':
    unittest.main()