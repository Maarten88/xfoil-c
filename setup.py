# -*- coding: utf-8 -*-
#   Copyright (c) 2019 D. de Vries
#
#   This file is part of XFoil, licensed under the GNU General Public License v3 or later.
#   See LICENSE for details.
"""Build script for xfoil_c.

CMake builds the plain C shared library libxfoil_c (not a CPython extension module)
and places it inside the xfoil_c package directory, where xfoil_c/_lib.py loads it
with ctypes. This mirrors xfoil-python/setup.py, but for Linux only for now.

    pip install .                        # build and install
    python setup.py build_ext --inplace  # build into ./xfoil_c/ for in-tree tests
    python setup.py build_ext --inplace --debug

Build variants (full description: header of CMakeLists.txt). Two switches, given as environment
variables or as -D options; setup.py always passes both to CMake, so a cached build directory never
keeps an old choice:
    XFOIL_ARCH   ""      (default) portable, bit-identical to the Fortran reference
                 native  this machine's CPU incl. FMA (~10-15 % faster, last-digit differences)
                 <cpu>   any gcc -march value, e.g. x86-64-v3
    XFOIL_FAST   OFF     (default) exact code paths
                 ON      source-level fast paths (#if XFOIL_FAST), e.g. the vectorised blsolv

    python setup.py build_ext --inplace                                      # exact (default)
    XFOIL_ARCH=native XFOIL_FAST=ON python setup.py build_ext --inplace      # fastest on this machine
    python setup.py build_ext --inplace -DXFOIL_ARCH=x86-64-v3 -DXFOIL_FAST=ON
    python -c "import xfoil_c; print(xfoil_c.build_info())"                  # which variant is loaded
"""
import glob
import os
import re
import shutil
import subprocess
import sys

from setuptools import setup
from setuptools.command.build_ext import build_ext
from setuptools.extension import Extension

__version__ = re.findall(
    r"""__version__ = ["']+([0-9\.]*)["']+""",
    open('xfoil_c/__init__.py').read(),
)[0]

LIB_PATTERNS = ('libxfoil_c.so', 'libxfoil_c.dylib', 'xfoil_c.dll', 'libxfoil_c.dll')

debug = '--debug' in sys.argv
if debug:
    sys.argv.remove('--debug')

# Command line -D flags are forwarded to CMake
cmake_cmd_args = [f for f in sys.argv if f.startswith('-D')]
for f in cmake_cmd_args:
    sys.argv.remove(f)


class CMakeExtension(Extension):

    def __init__(self, name, cmake_list_dir='.', **kwargs):
        super().__init__(name, sources=[], **kwargs)
        self.cmake_lists_dir = os.path.abspath(cmake_list_dir)


class CMakeBuild(build_ext):

    def build_extensions(self):
        try:
            subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError('Cannot find CMake executable')

        for ext in self.extensions:
            # Directory of the package the library must land in
            extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
            cfg = 'Debug' if debug else 'Release'

            cmake_args = [
                '-DCMAKE_BUILD_TYPE=%s' % cfg,
                '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir),
                '-DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), self.build_temp),
            ]
            for var, default in (('XFOIL_ARCH', ''), ('XFOIL_FAST', 'OFF')):
                if not any(a.startswith('-D%s=' % var) for a in cmake_cmd_args):
                    cmake_args.append('-D%s=%s' % (var, os.environ.get(var, default)))
            cmake_args += cmake_cmd_args

            os.makedirs(self.build_temp, exist_ok=True)
            subprocess.check_call(['cmake', ext.cmake_lists_dir] + cmake_args, cwd=self.build_temp)
            subprocess.check_call(['cmake', '--build', '.', '--config', cfg], cwd=self.build_temp)

    def copy_extensions_to_source(self):
        # --inplace: setuptools would copy a CPython extension named after ext.name, which
        # does not exist. Copy the plain shared library CMake built into the source package.
        build_py = self.get_finalized_command('build_py')
        for ext in self.extensions:
            package = '.'.join(ext.name.split('.')[:-1])
            src_dir = build_py.get_package_dir(package)
            built_dir = os.path.join(self.build_lib, *package.split('.'))
            libs = [f for pattern in LIB_PATTERNS for f in glob.glob(os.path.join(built_dir, pattern))]
            if not libs:
                raise RuntimeError('CMake build produced no library in ' + built_dir)
            for lib in libs:
                shutil.copy2(lib, src_dir)


def readme():
    with open('README.md') as f:
        return f.read()


setup(
    name='xfoil_c',
    version=__version__,
    description='C port of the stripped-down XFOIL python module',
    long_description=readme(),
    long_description_content_type='text/markdown',
    url='https://github.com/Maarten88/xfoil-c',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: C',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering',
    ],
    keywords='xfoil airfoil aerodynamic analysis',
    license='GPL-3.0-or-later',
    packages=['xfoil_c'],
    ext_modules=[CMakeExtension('xfoil_c._libxfoil_c')],
    cmdclass={'build_ext': CMakeBuild},
    install_requires=['numpy'],
    zip_safe=False,
)
