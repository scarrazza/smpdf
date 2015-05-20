# -*- coding: utf-8 -*-
"""
Created on Sun May 10 01:03:24 2015

@author: zah
"""
import unittest

from smpdflib.core import PDF
from smpdflib import config


class TestConfig(unittest.TestCase):

    def _test_bad_config(self, s):
        self.assertRaises(config.ConfigError, config.Config.from_yaml,s)

    def test_bad_spacing(self):
        s = (
"""observables:
   - {name: data/applgrid/atlas-incljets-r06-arxiv-1112.6297-eta7.root, order: 1 }
   - {name: data/applgrid/ttbar-xsectot-8tev.root, order: 1}
fmt: svg
actions:
   -savedata"""
         )
        self._test_bad_config(s)

    def test_default_propagation(self):
        s=  (
"""observables:
   - {name: data/applgrid/atlas-incljets-r06-arxiv-1112.6297-eta7.root, order: 1 }
   - {name: data/applgrid/ttbar-xsectot-8tev.root, order: 1}
fmt: svg
base_pdf: MMHT2014nlo68cl_1000MC
actions:
   - savedata

actiongroups:
   - pdfsets :
      - NNPDF30*_as_????
      - NNPDF30*_as_????_nf_?
      - NNPDF30_nlo_as_0119_atlas
      - MMHT2014nlo68cl_1000MC
     actions:
      - testas
      - asplots
      - savedata
     prefix: asgroup
   - pdfsets :
      - NNPDF30_nlo_as_0118
      - {name: MMHT2014nlo68cl_1000MC}
     base_pdf : NNPDF30_nlo_as_0118
     fmt: png
     actions:
      - violinplots
      - savedata
   - prefix: empty_action
     pdfsets :
      - NNPDF30*_as_????
      - NNPDF30*_as_????_nf_?
      - MMHT2014nlo68cl_1000MC"""
        )
        conf = config.Config.from_yaml(s)
        self.assertEqual(conf.actiongroups[2]['actions'], set(['exporthtml',
                                                               'exportcsv']))
        self.assertEqual(conf.actiongroups[0]['base_pdf'],
                         PDF('MMHT2014nlo68cl_1000MC'))

    def test_bad_name(self):
        s=  (
"""observables:
   - {name: data/applgrid/atlas-incljets-r06-arxiv-1112.6297-eta7.root, order: 1 }
   - {name: data/applgrid/ttbar-xsectot-8tev.root, order: 1}
fmt: svg
actions:
   - savedata

actiongroups:
   - pdfsets :
      - BADSETblablabla
      - NNPDF30*_as_????_nf_?
      - NNPDF30_nlo_as_0119_atlas
     actions:
      - testas
      - asplots
      - savedata
     prefix: asgroup
   - pdfsets :
      - NNPDF30_nlo_as_0118
      - {name: MMHT2014nlo68cl_1000MC}
     base_pdf : NNPDF30_nlo_as_0118
     fmt: png
     actions:
      - violinplots
      - savedata
   - prefix: empty_action
     pdfsets :
      - NNPDF30*_as_????
      - NNPDF30*_as_????_nf_?"""
        )
        self._test_bad_config(s)
        s=  (
"""observables:
   - {name: data/applgrid/atlas-incljets-r06-arxiv-1112.6297-eta7.root, order: 1 }
   - {name: data/applgrid/ttbar-xsectot-8tev.root, order: 1}
fmt: svg
actions:
   - savedata

actiongroups:
   - pdfsets :
      - BADSETblablabla
     actions:
      - testas
      - asplots
      - savedata
     prefix: asgroup
   - pdfsets :
      - NNPDF30_nlo_as_0118
      - {name: MMHT2014nlo68cl_1000MC}
     base_pdf : NNPDF30_nlo_as_0118
     fmt: png
     actions:
      - violinplots
      - savedata
   - prefix: empty_action
     pdfsets :
      - NNPDF30*_as_????
      - NNPDF30*_as_????_nf_?"""
        )
        self._test_bad_config(s)
        s=  (
"""observables:
   - {name: data/applgrid/atlas-incljets-r06-arxiv-1112.6297-eta7.root, order: 1 }
   - {name: data/applgrid/ttbar-xsectot-8tev.root, order: 1}
fmt: svg
actions:
   - savedata

actiongroups:
   - pdfsets :
      - BADSET??blablabla*
      - NNPDF30*_as_????_nf_?
      - NNPDF30_nlo_as_0119_atlas
     actions:
      - testas
      - asplots
      - savedata
     prefix: asgroup
   - pdfsets :
      - NNPDF30_nlo_as_0118
      - {name: MMHT2014nlo68cl_1000MC}
     base_pdf : NNPDF30_nlo_as_0118
     fmt: png
     actions:
      - violinplots
      - savedata
   - prefix: empty_action
     pdfsets :
      - NNPDF30*_as_????
      - NNPDF30*_as_????_nf_?"""
        )
        self._test_bad_config(s)
        s=  (
"""observables:
   - {name: data/applgrid/atlas-incljets-r06-arxiv-1112.6297-eta7.root, order: 1 }
   - {name: data/applgrid/ttbar-xsectot-8tev.root, order: 1}
fmt: svg
actions:
\t- savedata
"""
        )
        self._test_bad_config(s)


if __name__ == '__main__':
    unittest.main()