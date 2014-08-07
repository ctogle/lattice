from math import log
import numpy as np
import itertools as it

import random as rm

import pdb

def set_state(state, time, latt):
	state[0] = time
	state[1] = latt.total_population()

def capture(data, state, capture_dex, target_dexes):
	data[:,capture_dex] = [state[dex] for dex in target_dexes]

class lattice_site(object):

	def __init__(self, pos, pop = 3, prop = 0.0):
		self.position = pos
		self.population = pop
		self.prop = prop

	def propensity(self):
		return self.prop

	def update_propensity(self, latt):
		self.prop = 1.0

	def scale_propensity(self, scale):
		self.prop *= scale

	def react(self):
		#print 'site at', self.position, 'reacted!'
		self.population += 1

class lattice(object):

	def __init__(self, *args, **kwargs):
		#self.kwargs = kwargs
		#self.args = args
		if 'dims' in kwargs.keys(): dims = kwargs['dims']
		else: dims = 2
		if 'ax_length' in kwargs.keys(): ax_length = kwargs['ax_length']
		else: ax_length = 10
		self.make_lattice(dims, ax_length)

	def make_lattice(self, dims, maxax):
		self.dimensions = dims
		self.max_axis_size = maxax
		v_site = np.vectorize(lattice_site)
		init_arry = np.empty((maxax,)*dims,dtype = object)
		self.locs = [x for x in it.product(*[range(maxax)]*dims)]
		for pos in self.locs: init_arry[pos] = pos
		lattice = np.empty((maxax,)*dims,dtype = object)
		lattice[:] = v_site(init_arry)
		self._lattice_ = lattice

	def total_population(self):
		pop = 0
		for pos in self.locs: pop += self._lattice_[pos].population
		return pop

	def update_propensities(self):
		for pos in self.locs:
			self._lattice_[pos].update_propensity(self._lattice_)

	def react(self, ldex):
		if ldex == -1: self.update_propensities()
		else:
			self._lattice_[self.locs[ldex]].react()
			self.update_propensities()

	def total_propensity(self):
		prop = 0
		for pos in self.locs: prop += self._lattice_[pos].propensity()
		return prop

	def scale_propensities(self, scl):
		for pos in self.locs: self._lattice_[pos].scale_propensity(scl)

	def pick_lattice_site(self):
		ax_rand = rm.randrange(len(self.locs))
		return ax_rand

	def population_grid(self):
		def fix(d):
			if d in relev: return slice(None)
			else: return fixins[d]
		fixins = [0 for d in range(self.dimensions)]
		relev = [0,1]
		fixed = [fix(d) for d in range(len(fixins))]
		maxax = self.max_axis_size
		surf = np.zeros((maxax, maxax), dtype = np.int)
		surf[:,:] = [[sit.population for sit in ax] 
					for ax in self._lattice_[fixed]]
		return surf

def_string = ''
def simulate(sys_string = def_string):
	_lattice_ = lattice(dims = 2, ax_length = 100)

	time = 0.0
	last_time = 0.0
	end_time = 1.0
	incr_time = 0.02

	total_captures = end_time/incr_time
	capture_count = 0

	targets = ['time', 'population']
	plot_targets = ['time', 'population']
	target_dexes = [targets.index(ta) for ta in plot_targets]

	maxax = _lattice_.max_axis_size
	state = np.zeros(shape = (len(targets)), dtype = np.float)
	data = np.zeros(shape = (len(plot_targets), total_captures), dtype = np.float)
	pop_surfaces = np.zeros(shape = (total_captures, maxax, maxax), dtype = np.float)
	set_state(state, time, _lattice_)

	iteration = 0
	while capture_count < total_captures:
		iteration += 1
		propensity_total = _lattice_.total_propensity()

		if propensity_total > 0.0:
			propensity_total_inv = 1.0/propensity_total
			time_step = -1.0 * log(rm.random()) * propensity_total_inv
			_lattice_.scale_propensities(propensity_total_inv)
			lsite_dex = _lattice_.pick_lattice_site()
		else:
			time_step = incr_time
			lsite_dex = -1

		time += time_step
		set_state(state, time, _lattice_)

		real_time = state[0]
		while last_time < real_time and capture_count < total_captures:
			state[0] = last_time
			last_time += incr_time
			capture(data, state, capture_count, target_dexes)
			pop_surfaces[capture_count] = _lattice_.population_grid()
			capture_count += 1
			print 'time', last_time
		state[0] = real_time

		_lattice_.react(lsite_dex)

	return data,pop_surfaces,targets+['population_surfaces']

if __name__ == '__main__':
	data, targets, pop_surfaces, iteration = simulate()
	pdb.set_trace()

