import unittest

import modular_core.libfundamental as lfu
from modular_core.libsimcomponents import ensemble_manager as mng

import os, sys, pdb

#log = open(os.path.join(os.getcwd(), 'test_ensemble.log'), 'w')
#sys.stdout = log

import lattice as mod

class dummyTestCase(unittest.TestCase):
	"""Tests for `dummy module`."""
	simple_mcfg = os.path.join(os.getcwd(), 
				'lattice_dep_mcfgs', 
				'lattice_example.mcfg')
	mn = mng()
	en = mn.add_ensemble(module = mod.main.module_name)

	def pause(self, *args, **kwargs):
		sys.stdout = sys.__stdout__
		pdb.set_trace()
		sys.stdout = log

	def test_can_make_ensemble(self):
		"""module successfully imported?"""
		self.assertFalse(mod.main == None)
		mod_name = mod.main.module_name
		ensem = self.mn.add_ensemble(module = mod_name)
		self.assertFalse(ensem is None)

	def test_can_run_mcfg(self):
		ran = self.en.run_mcfg(self.simple_mcfg)
		out = self.en.produce_output()
		self.assertTrue(ran)
		self.assertTrue(out)

if __name__ == '__main__':
    unittest.main()

