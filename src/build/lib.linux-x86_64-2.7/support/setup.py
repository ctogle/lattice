from subprocess import call

import Cython.Compiler.Options as copts
copts.annotate = True

from distutils.core import setup
#from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy

ext_modules = [Extension('lattice_extension', ['lattice_extension.pyx'])]
setup(
	name = 'lattice_util', 
	ext_modules = ext_modules, 
	cmdclass = {'build_ext': build_ext}, 
	include_dirs = [numpy.get_include()]
)

call(['cp', 'liblattice_sim.py', 'liblattice_sim.pyx'])

ext_modules = [Extension('lattice_simulator', ['liblattice_sim.pyx'])]
setup(
	name = 'lattice_simulator', 
	ext_modules = ext_modules, 
	cmdclass = {'build_ext': build_ext}, 
	include_dirs = [numpy.get_include()]
)

#python setup.py build_ext --inplace


