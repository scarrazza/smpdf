# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
Created on Fri Jan 23 12:11:23 2015

@author: zah
"""
from __future__ import print_function
import os
import os.path as osp
import re
import subprocess
import fnmatch

import yaml


_indexes_to_names = None
_names_to_indexes = None
#TODO Use Python 3 cache when sensible
infofiles = {}

def expand_names(globstr):
    return fnmatch.filter(get_names_to_indexes().keys(), globstr)
        
def get_indexes_to_names():
    global _indexes_to_names
    if _indexes_to_names is None:
        _indexes_to_names = parse_index(get_index_path())
    return _indexes_to_names

def get_names_to_indexes():
    global _names_to_indexes
    if _names_to_indexes is None:
        _names_to_indexes = {name:index for index,name in 
                            get_indexes_to_names().items()}
    return _names_to_indexes

def get_pdf_indexes(name):
    """Get index in the amc@nlo format"""
    info = parse_info(name)
    ind = info['SetIndex']
    num_members = info['NumMembers']
    return {'lhapdf_id' : ind, 
            'lhapdf_min': ind + (num_members > 1), 
            'lhapdf_max': ind + num_members - 1}

def get_pdf_name(index):
    return get_indexes_to_names()[str(index)]

def parse_index(index_file):
    d = {}
    name_re = '(\d+)\s+(\S+)'
    with open(index_file) as localfile:
        for line in localfile.readlines():
            m = re.match(name_re, line)
            index = m.group(1)
            d[index] = m.group(2)
    return d

def parse_info(name):
    global infofiles
    if name in infofiles: 
        return infofiles[name]
    infofilename = osp.join(get_lha_path(), name, name + '.info')
    with open(infofilename) as infofile:
        result = yaml.load(infofile)
    infofiles[name] = result
    return result

def get_lha_path():
    lhapdf_folder = subprocess.check_output(['lhapdf-config', '--datadir'])
    lhapdf_folder = lhapdf_folder.decode().rstrip()
    return lhapdf_folder

def get_index_path(folder = None):
    if folder:
        index_file = os.path.join(folder, 'pdfsets.index')
    if folder is None or not osp.exists(index_file):
        folder = get_lha_path()
    index_file = os.path.join(folder, 'pdfsets.index')
    return index_file
