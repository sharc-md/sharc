#******************************************
#
#    SHARC Program Suite
#
#    Copyright (c) 2023 University of Vienna
#
#    This file is part of SHARC.
#
#    SHARC is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    SHARC is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    inside the SHARC manual.  If not, see <http://www.gnu.org/licenses/>.
#
#******************************************


import os
import sys


# py_setuptools applies some customization to the setuptools code
# Therefore, we load all stuff from setuptools through py_setuptools
from setuptools import setup, Extension#, settings

# 
#from pysharc_setup import pysharc_extension











# path to the pysharc files
# (relative to the directory from where the setup script is started)
pysharc_path = 'pysharc_src'

# Files with C code
# inside pysharc_path/
pysharc_cfiles = ['pysharc.c', 'pysharc_tools.c']



# define Libraries
mkl_libs = []
basic_libs = ['sharc']#, 'hdf5', 'hdf5_hl', 'netcdf']
#basic_libs = ['sharc', 'gfortran', 'hdf5', 'hdf5_hl', 'netcdf']
extra_compile_args = ['-std=c99', '-Wall',]
#extra_compile_args += ['-D__PYTHON_DEBUG__', '-Wall']






pysharc_extension = Extension('sharc/sharc',
            [ os.path.join(pysharc_path, fname) for fname in pysharc_cfiles ],
            include_dirs = [ 'include', ],
            library_dirs = [ 'lib', os.environ.get('ANACONDA', '') + '/lib' ],
            libraries = mkl_libs + basic_libs,
            extra_compile_args = extra_compile_args,
)







prog_name = 'sharc'
version = '0.1.1.'
description='Python API for the SHARC MD Package'
author = 'Maximilian F.S.J. Menger'
author_email='maximilian.menger@univie.ac.at'
url='https://sharc-md.org'


sharc_packages = [ 'sharc', 'sharc.pysharc' ]
sharc_scripts =  [ 'bin/pysharc_lvc.py', 'bin/pysharc_qmout.py' ]


setup(name=prog_name,
#      version=version,
      description=description,
      author=author,
      author_email=author_email,
      url=url,
      packages = sharc_packages,
      scripts = sharc_scripts,
      ext_modules=[pysharc_extension],
      install_requires=[
#                    'netcdf4', # for correct libraries
#                    'hdf5',    # for the libraries
                    'mkl',
                          ],
)
