# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 15:41:26 2015

@author: zah
"""
import sys
import os
from contextlib import contextmanager, redirect_stdout

@contextmanager
def supress_stdout():
    devnull = open(os.devnull, 'wb')
    try:
        stdout_flieno = sys.stdout.fileno()
    except ValueError:
        redirect = False
    else:
        redirect = True
        sys.stdout.flush()
        #sys.stdout.close()
        devnull_fileno = devnull.fileno()
        saved_stdout_fd = os.dup(stdout_flieno)
        os.dup2(devnull_fileno, stdout_flieno)

    with redirect_stdout(devnull):
        yield
    if redirect:
        os.dup2(stdout_flieno, saved_stdout_fd)