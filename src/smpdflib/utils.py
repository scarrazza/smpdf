# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 15:41:26 2015

@author: zah
"""
import sys
import os
from contextlib import contextmanager, redirect_stdout

@contextmanager
def _supress_stdout():
    """This fits the 90% use case. It will fail badly inside ipython
    or when sys.stdout is a strange object."""
    sys.stdout.flush() # <--- important when redirecting to files
    # Duplicate stdout (file descriptor 1)
    # to a different file descriptor number
    # /dev/null is used just to discard what is being printed
    devnull = os.open('/dev/null', os.O_WRONLY)
    newstdout = os.dup(1)
    # Duplicate the file descriptor for /dev/null
    # and overwrite the value for stdout (file descriptor 1)
    os.dup2(devnull, 1)
    try:

        yield

    finally:
        # Close devnull after duplication (no longer needed)
        os.close(devnull)

        # Use the original stdout to still be able
        # to print to stdout within python
        sys.stdout = os.fdopen(newstdout, 'w')

def supress_stdout():
    try:
        fd = sys.stdout.fileno()
    except ValueError:
         fd = None
    if fd == 1:
        return _supress_stdout()
    else:
        return redirect_stdout(open(os.devnull, 'w'))