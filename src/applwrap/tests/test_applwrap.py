# -*- coding: utf-8 -*-
"""
Created on Fri May 15 16:45:30 2015

@author: zah
"""
import unittest

import applwrap


class TestConfig(unittest.TestCase):
    def test_bad_values(self):
        with self.assertRaises(ValueError):
            applwrap.initobs("patata")
        with self.assertRaises(ValueError):
            applwrap.initpdf("patata")

if __name__ == '__main__':
    unittest.main()
