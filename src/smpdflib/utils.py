# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 15:41:26 2015

@author: zah
"""
import sys
import os
import multiprocessing
import logging
import atexit
from logging.handlers import QueueHandler, QueueListener
from contextlib import contextmanager, redirect_stdout

import pandas as pd

import applwrap

@contextmanager
def _supress_stdout():
    """This fits the 90% use case. It will fail badly inside ipython
    or when sys.stdout is a strange object."""
    sys.stdout.flush() # <--- important when redirecting to files
    # Duplicate stdout (file descriptor 1)
    # to a different file descriptor number
    # /dev/null is used just to discard what is being printed
    devnull = os.open(os.devnull, os.O_WRONLY)
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

#https://gist.github.com/Zaharid/d4f8c9a44ce7941b0c37
__queue = None
def get_logging_queue():
    #This probably will have to be refactored to get access to manager as well.
    global __queue, __manager
    if __queue is None:
        m = multiprocessing.Manager()
        __manager = m
        q = m.Queue(-1)
        #https://docs.python.org/3/howto/logging-cookbook.html
        listener = QueueListener(q, *logging.getLogger().handlers)
        listener.start()
        def exithandler():
            q.join()
            listener.stop()
            #Seems to help silencing bugs...
            import time; time.sleep(0.2)
            m.shutdown()

        atexit.register(exithandler)
        __queue = q
        return q
    return __queue

def initlogging(q, loglevel):
    handler = QueueHandler(q)
    l = logging.getLogger()
    l.level = loglevel
    l.handlers = [handler,]
    if not l.isEnabledFor(logging.DEBUG):
        applwrap.setverbosity(0)
    atexit.register(lambda: q.join())

def save_html(df, path):
    import jinja2

    env = jinja2.Environment(loader = jinja2.PackageLoader('smpdflib',
                                                           'templates'))
    template = env.get_template('results.html.template')
    def remark_formatter(remarks):
        if not remarks:
            return ''
        else:
            return '<ul>%s</ul>' % '\n'.join('<li>%s</li>' %
                   jinja2.escape(remark) for remark in remarks)

    #http://stackoverflow.com/questions/26277757/pandas-to-html-truncates-string-contents
    with pd.option_context('display.max_colwidth', -1):
        table = df.to_html(
                             formatters={'Remarks':remark_formatter},
                             escape = False)
    result = template.render(table=table)
    with open(path, 'w') as f:
        f.write(result)