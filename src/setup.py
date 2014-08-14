#!/usr/bin/python
from os.path import isfile, join
import glob
import os
import re

from setuptools import setup, Extension

if isfile("MANIFEST"):
    os.unlink("MANIFEST")

TOPDIR = os.path.dirname(__file__) or "."
#VERSION = re.search('__version__ = "([^"]+)"',
#	open(TOPDIR + "/src/pkg1/__init__.py").read()).group(1)
VERSION = '1.0'

core_modules = [
	'lattice.scripts.liblattice', 

	'support.liblattice_sim', 
				]

#def ignore_res(f):
#	#if f.startswith('__') or f.startswith('_.'): return True
#	#else: return False
#	return False
#res_dir = 'lattice/mcfgs/'
#res_fis = [f for f in os.listdir(os.path.join(
#	os.getcwd(), 'lattice', 'mcfgs')) if not ignore_res(f)]
#res_files = [res_dir + f for f in res_fis]

#requirements = [
#	'modular_core >= 1.0', 
#			]

setup(
	name="lattice-module-pkg",
	version = VERSION,
	description = "lattice-module pkg",
	author = "ctogle",
	author_email = "cogle@vt.edu",
	url = "http://github.com/ctogle/lattice",
	license = "MIT License",
	long_description =
"""\
This is the lattice module package of modular
""",
	#install_requires = requirements, 
	packages = ['lattice', 'lattice.scripts', 'support'], 
	py_modules = core_modules, 
	#dependency_links = [], 
	zip_safe = False,
	#ext_package = 'modular_core.modules.chemicallite_support',
	#data_files=[('lattice/mcfgs', res_files)], 
	ext_modules = [],
	)




