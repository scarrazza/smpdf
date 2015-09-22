from __future__ import print_function
import sys
from setuptools import setup, Extension, find_packages
import subprocess
import platform

if sys.version_info < (3,4):
    print("SMPDF requires Python 3.4 or later", file=sys.stderr)
    sys.exit(1)

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
extra_compile_args = lhapdf_includes + applgrid_includes
extra_link_args = lhapdf_includes + applgrid_includes + applgrid_libs + lhapdf_libs

if platform.system() == 'Darwin':
    mac_ver = platform.mac_ver()[0]
    extra_compile_args += ['-mmacosx-version-min=%s' % mac_ver]
    extra_link_args += ['-mmacosx-version-min=%s' % mac_ver]
module1 = Extension('applwrap',
                    extra_compile_args = extra_compile_args ,
                    extra_link_args = extra_link_args,
                    sources = ['src/applwrap/applwrap.cc'], language="c++")

setup (name = 'smpdf',
       version = '1.0rc1',
       description = "PDFs for phenomenology",
       author = 'Stefano Carrazza and Zahari Kassabov',
       author_email = 'kassabov@to.infn.it',
       url = 'https://github.com/scarrazza/smpdf',
       long_description = "See `smpdf --help` for the full list of options",
       ext_modules = [module1],
       #http://stackoverflow.com/questions/24799146/use-multiprocessing-process-inside-a-script-installed-by-setuptools
       entry_points = {'console_scripts': ['smpdf = _smpdf_scripts.smpdf:main']},
       package_dir = {'': 'src'},
       packages = find_packages('src'),
       package_data = {
            '':['*.template', '*.mplstyle']
       },
       zip_safe = False,
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
