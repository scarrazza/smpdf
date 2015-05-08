from setuptools import setup, Extension
import subprocess

def call_command(command):
    l = command.split()
    result = subprocess.check_output(l)
    return result.decode().rstrip()

lhapdf_includes = call_command('lhapdf-config --cflags').split()
lhapdf_libs = call_command('lhapdf-config --libs').split()

applgrid_includes = call_command('applgrid-config --cxxflags').split()
applgrid_libs = call_command('applgrid-config --ldflags').split()
module1 = Extension('applwrap',
                    extra_compile_args = (lhapdf_includes + applgrid_includes +
                                          applgrid_libs + lhapdf_libs),
                    extra_link_args = (lhapdf_includes + applgrid_includes  +
                                          applgrid_libs + lhapdf_libs),
                    sources = ['applwrap/applwrap.cc'], language="c++")

setup (name = 'PackageName',
       version = '1.0',
       description = 'This is a demo package',
       author = 'Martin v. Loewis',
       author_email = 'martin@v.loewis.de',
       url = 'https://docs.python.org/extending/building',
       long_description = '''
This is really just a demo package.
''',
       ext_modules = [module1],
       scripts = ['smpdf.py']
       )
