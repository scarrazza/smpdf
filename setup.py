from __future__ import print_function
import sys
from setuptools import setup, Extension, find_packages
import subprocess

def call_command(command):
    l = command.split()
    try:
        result = subprocess.check_output(l)
    except OSError as e:
        print("Could not call %s: %s.\n"
              "Please make sure the relevant command is isntalled."
              % (l[0], e.strerror), file=sys.stderr)
        sys.exit(1)
    return result.decode().rstrip()

lhapdf_includes = call_command('lhapdf-config --cflags').split()
lhapdf_libs = call_command('lhapdf-config --libs').split()

applgrid_includes = call_command('applgrid-config --cxxflags').split()
applgrid_libs = call_command('applgrid-config --ldflags').split()
module1 = Extension('applwrap',
                    extra_compile_args = (lhapdf_includes + applgrid_includes
                                           ),
                    extra_link_args = (lhapdf_includes + applgrid_includes  +
                                          applgrid_libs + lhapdf_libs),
                    sources = ['src/applwrap/applwrap.cc'], language="c++")

setup (name = 'smpdf',
       version = '0.1',
       description = "PDFs for phenomenology",
       author = 'Stefano Carrazza and Zahari Kassabov',
       author_email = 'kassabov@to.infn.it',
       url = 'https://github.com/scarrazza/smpdf',
       long_description = "See `smpdf --help` for the full list of options",
       ext_modules = [module1],
       scripts = ['scripts/smpdf.py'],
       package_dir = {'': 'src'},
       packages = find_packages('src'),
       classifiers=[
            'License :: OSI Approved :: BSD License',
            'Operating System :: OS Independent',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Physics',
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.4',
            ],
       )
