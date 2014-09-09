#cython: profile=True,boundscheck=False

from lattice_extension import validate_rxn as validate_rxn
from lattice_extension import propensity as rxn_propensity
from lattice_extension import update_propensities as agents_propensity
from lattice_extension import update_agents_react_prop
from lattice_extension import choose_biased_index as cbindex
from lattice_extension import pick_agent_to_act as pick_agent

from math import log
import numpy as np
import itertools as it
import time as ti

import random as rm

import modular_core.libprofile as lprf

import pdb

def set_state(state, surface_state, time, latt):
    state[0] = time
    state[1] = latt.agent_count
    state[2] = np.mean(surface_state[0])
    #state[1] = latt.total_population()
    #state[2:] = spectrum_population(
    #    latt.spec_count, latt.agent_count, 
    #    latt.species_index_dict, latt.agents)
    state[3:] = latt.spectrum_population()

def set_surface_state(surf_state, time, latt):
    surf_state[0] = latt.population_grid()
    surf_state[1] = latt.identity_grid()
    if latt.has_resources:
        surf_state[2:] = latt.spectrum_resource_grid()

def capture(data, state, capture_dex, target_dexes):
    data[:,capture_dex] = [state[dex] for dex in target_dexes]

class agent(object):

    def __init__(self, variety, motion_prop = 0.0, react_prop = 0.0):
        self.variety = variety
        self.motion_prop = motion_prop
        self.react_prop = react_prop
        self.total_prop = self.react_prop + self.motion_prop
        self.behaviors = [self.move, self.react]

    def move_to(self, pos, lattice):
        site = lattice._lattice_[self.pos]
        site.agents.remove(self)
        newpopulace = site.populace.replace(self.variety, '', 1)
        site.populace = newpopulace
        site.population = len(newpopulace)
        site.who = site.identity()
        lattice.flagged_locs.append(self.pos)

        self.pos = pos
        newsite = lattice._lattice_[pos]
        newsite.agents.append(self)
        newsite.populace += self.variety
        newsite.population = len(newsite.populace)
        newsite.who = newsite.identity()
        lattice.flagged_locs.append(pos)

    def act(self, lattice):
        prtot = self.total_prop
        if prtot > 0.0:
            #normed = [self.motion_prop/prtot, self.react_prop/prtot]
            #behvs = [normed[0], normed[0] + normed[1]]
            #rd = rm.random()
            #if behvs[0] > rd:
            #    self.move(lattice)
            #elif behvs[1] > rd:
            #    self.react(lattice)
            #else: pdb.set_trace()
            behvs = [self.motion_prop, self.react_prop]
            behv = cbindex(behvs,2)
            self.behaviors[behv](lattice)
        else: print 'PASSING!'

    def move(self, lattice):
        new_pos = lattice._lattice_[self.pos].diffuse_pressure()
        #new_pos = lattice._lattice_[self.pos].diffuse_random()
        self.move_to(new_pos, lattice)

    def react(self, lattice):
        site = lattice._lattice_[self.pos]
        populace = site.populace
        resources = site.resources
        valid = [rx for rx in lattice.reactions	
            if validate_rxn(rx.used, populace, resources)]
        rxn_props = [rxn_propensity(
            rx.used, populace, rx.rate) 
                for rx in valid]
        rxdex = pick_rxn(rxn_props)
        rxn = valid[rxdex]
        rxn.react(site, lattice)
        lattice.flagged_locs.append(self.pos)

    def update_propensity(self, lattice):
        lattice.total_agent_propensity -= self.total_prop
        self.update_react_prop(lattice)
        self.update_motion_prop(lattice)
        self.total_prop = self.react_prop + self.motion_prop
        lattice.total_agent_propensity += self.total_prop

    def update_react_prop(self, lattice):
        site = lattice._lattice_[self.pos]
        rxns = lattice.reactions
        self.react_prop = update_agents_react_prop(
            site, rxns, lattice.reaction_count)

    def update_motion_prop(self, lattice):            
        self.motion_prop = lattice._lattice_[
            self.pos].population**lattice.dimensions

class agent_immobile(agent):

    def move_to(self, pos, lattice):
        print 'immobile should not be asked to move!'

    def act(self, lattice):
        self.react(lattice)

    def update_propensity(self, lattice):
        lattice.total_agent_propensity -= self.total_prop
        self.update_react_prop(lattice)
        self.total_prop = self.react_prop
        lattice.total_agent_propensity += self.total_prop

    def update_motion_prop(self, lattice):            
        print 'immobile agent should not be asked up update motion_prop!!'
        self.motion_prop = 0.0

def pick_rxn(rxn_props):
    table = []
    ptot = sum(rxn_props)
    rd = rm.uniform(0.0, ptot)
    for p in rxn_props: table.append(sum(table) + p)
    for tx, t in enumerate(table):
        if rd < t: return tx

class reaction(object):

    def __init__(self, rxn, specs):
        def convert(ag):
            stoch = int(ag[ag.find('(')+1:ag.find(')')])
            spec = ag[ag.find(')')+1:]
            isresource = not spec in specs
            return (stoch, spec, isresource)
        self.rxnstr = rxn
        rxnspl = rxn.split('->')
        front = rxnspl[0]
        self.rate = float(front[front.find('[')+1:front.find(']')])
        used, prod = front[:front.find('[')], rxnspl[1]
        if not used == '': used = used.split('+')
        else: used = []
        if not prod == '': prod = prod.split('+')
        else: prod = []
        self.used = [convert(us) for us in used]
        self.prod = [convert(pr) for pr in prod]
        self.flux = len(self.prod) - len(self.used)

    def react(self, site, lattice):
        pop = site.populace[:]
        resources = site.resources
        res_varis = [res.variety for res in resources]
        for u in self.used:
            if u[2]:
                rdex = res_varis.index(u[1])
                resource = resources[rdex]
                resource.quantity -= u[0]
            else: pop = pop.replace(u[1], '', u[0])
        for p in self.prod:
            if p[2]:
                rdex = res_varis.index(p[1])
                resource = resources[rdex]
                resource.quantity += p[0]
            else: pop += p[1]*p[0]
        site.resolve_agents(lattice, pop)

class site_resource(object):
    def __init__(self, pos, quantity = 1, variety = 'food'):
        self.variety = variety
        self.quantity = quantity

#site_topology instances describe the connectivity of an N
# dimensional lattice site to its neighboring sites
class site_topology(object):
    def __init__(self, pos, lshape, max_occupancy = None):
        self.pos = pos
        self.lshape = lshape
        self.dims = len(pos)
        self.max_occupancy = max_occupancy
        self.set_step_positions()

    def set_step_positions(self):
        pos = self.pos
        self.dn_steps = [list(pos[:]) for x in range(self.dims)]
        self.up_steps = [list(pos[:]) for x in range(self.dims)]
        axlen = self.lshape[0]
        for d in range(self.dims):
            if pos[d] == 0:
                dstep = 1 - axlen
                ustep = 1
            elif pos[d] == axlen - 1:
                dstep = 1
                ustep = 1 - axlen
            else:
                dstep = 1
                ustep = 1

            self.dn_steps[d][d] -= dstep
            self.up_steps[d][d] += ustep

        self.up_steps = [tuple(x) for x in self.up_steps]
        self.dn_steps = [tuple(x) for x in self.dn_steps]
        self.all_steps = self.up_steps + self.dn_steps
        self.neighbor_count = len(self.all_steps)

    def step(self, axdex, dir_):
        if dir_: new = self.up_steps[axdex]
        else: new = self.dn_steps[axdex]
        return new

class lattice_site(object):

    def __init__(self, args, pop = ''):
        self.position = args[0]
        self.lshape = args[1]
        self.species = args[2].species
        self.resources = args[3]
        #self.resources = [site_resource(self.position)]
        self.resource_names = [r.variety for r in self.resources]
        self.dims = len(self.position)
        self.agents = []
        self.populace = pop
        self.population = len(pop)
        self.who = self.identity()
        self.lattice = args[2]
        max_occupancy = args[2].max_occupancy
        self.topology = site_topology(
            self.position, self.lshape, max_occupancy)
                
    def resolve_agents(self, lattice, pop):
        local_ags = [ag for ag in self.agents if not ag.variety == '_']
        #local_ags = [ag for ag in lattice.agents if ag.pos == self.position]
        local_pop = [ag.variety for ag in local_ags]
        left = []
        for ag in local_ags:
            if ag.variety in pop: pop = pop.replace(ag.variety, '', 1)
            else: left.append(ag)

        for p in pop: lattice.put_onto_lattice_at_site(self.position, agent(p))
        for p in left: lattice.remove_from_lattice_site(self.position, p)

    def diffuse_pressure(self):
        adjacents = self.topology.all_steps
        ncnt = self.topology.neighbor_count
        latt = self.lattice._lattice_
        adjac_pops = [latt[pos].population for pos in adjacents]
        max_po = max(adjac_pops)
        adjac_pops = [-1.0*po + max_po for po in adjac_pops]
        lookup = cbindex(adjac_pops, ncnt)
        step = adjacents[lookup]
        return step

    def diffuse_random(self):
        step_dim = rm.randrange(self.dims)
        dir_ = rm.random() > 0.5
        step = self.topology.step(step_dim, dir_)
        return step

    def identity(self):
        if not self.populace: return 0
        for sdx in range(len(self.species)):
            sp = self.species[sdx]
            if sp in self.populace: return sdx + 1
        return 0

    def what(self, resource_dex):
        res = self.resources[resource_dex]
        return res.quantity

class lattice(object):

    def __init__(self, *args, **kwargs):
        if 'dims' in kwargs.keys(): dims = kwargs['dims']
        else: dims = 2
        if 'ax_length' in kwargs.keys(): ax_length = kwargs['ax_length']
        else: ax_length = 5
        if 'species' in kwargs.keys(): specs = kwargs['species']
        else: specs = []
        if 'reactions' in kwargs.keys(): rxns = kwargs['reactions']
        else: rxns = []
        if 'birthing_flag' in kwargs.keys(): bflag = kwargs['birthing_flag']
        else: bflag = False
        if 'max_occupancy' in kwargs.keys():
            max_occupancy = kwargs['max_occupancy']
        else: max_occupancy = None
        if 'resources' in kwargs.keys(): ress = kwargs['resources']
        else: ress = []
        self.agents = []
        self.agent_count = 0
        self.species = specs
        self.spec_count = len(specs)
        self.total_agent_propensity = 0
        self.agents_propensities = []
        #spdxdict = {}
        #for sdx in xrange(self.spec_count):
        #    spec = self.species[sdx]
        #    spdxdict[spec] = sdx
        #self.species_index_dict = spdxdict
        self.species_spectrum = {}
        for ke in self.species: self.species_spectrum[ke] = 0
        self.reactions = [reaction(rx, specs) for rx in rxns]
        self.reaction_count = len(rxns)
        #self.birth_reactions =\
        #    [rx for rdx, rx in enumerate(self.reactions) if rdx in brxns]
        self.birth_flag = bflag
        self.max_occupancy = max_occupancy
        self.resources = ress
        if len(ress) == 0: self.has_resources = False
        else: self.has_resources = True
        self.make_lattice(dims, ax_length, ress)

    def make_lattice(self, dims, maxax, ress):
        #must add ress contents here; right now its hardcoded
        self.dimensions = dims
        self.max_axis_size = maxax
        v_site = np.vectorize(lattice_site)
        init_arry = np.empty((maxax,)*dims,dtype = object)
        self.locs = [x for x in it.product(*[range(maxax)]*dims)]
        lattice = np.empty((maxax,)*dims,dtype = object)
        shp = lattice.shape
        for pos in self.locs:
            res = []
            for re in ress:
                if re[0] == 'all' or re[0].count(pos) > 0:
                    res.append(site_resource(pos, re[2], re[1]))
            init_arry[pos] = (pos, shp, self, res)
        lattice[:] = v_site(init_arry)
        self._lattice_ = lattice
        self.flagged_locs = []
        if self.birth_flag:
            for pos in self.locs:
                bagent = agent_immobile('_')
                self.put_onto_lattice_at_site(pos, bagent)
                self.flagged_locs.append(pos)

    def put_onto_lattice_at_site(self, pos, agent):
        self.agents.append(agent)
        if pos == 'random':pos = rm.choice(self.locs)
        agent.pos = pos
        site = self._lattice_[pos]
        site.agents.append(agent)
        if not agent.variety == '_':
            self.species_spectrum[agent.variety] += 1
            site.populace += agent.variety
            self.agent_count += 1
        site.population = len(site.populace)
        site.who = site.identity()

    def remove_from_lattice_site(self, pos, agent):
        self.species_spectrum[agent.variety] -= 1
        site = self._lattice_[pos]
        site.populace = site.populace.replace(agent.variety, '', 1)
        site.population = len(site.populace)
        site.who = site.identity()
        site.agents.remove(agent)
        self.agents.remove(agent)
        self.agent_count -= 1

    def total_propensity(self):
        return self.total_agent_propensity

    def act_agent(self, agdex):
        agent = self.agents[agdex]
        agent.act(self)

    def get_clean_fixed_surf(self):
        def fix(d):
            if d in relev: return slice(None)
            else: return fixins[d]
        fixins = [0 for d in range(self.dimensions)]
        relev = [0,1]
        fixed = [fix(d) for d in range(len(fixins))]
        maxax = self.max_axis_size
        surf = np.zeros((maxax, maxax), dtype = np.int)
        return surf, fixed

    def population_grid(self):
        surf, fixed = self.get_clean_fixed_surf()
        surf[:,:] = [[sit.population for sit in ax] 
            for ax in self._lattice_[fixed]]
        return surf

    def identity_grid(self):
        surf, fixed = self.get_clean_fixed_surf()
        surf[:,:] = [[sit.who for sit in ax] 
            for ax in self._lattice_[fixed]]
        return surf

    def spectrum_population(self):
        spectrum = [0]*self.spec_count
        for sdx in xrange(self.spec_count):
            sp = self.species[sdx]
            spectrum[sdx] = self.species_spectrum[sp]
        return spectrum

    def spectrum_resource_grid(self):
        def resource_grid(res):
            surf, fixed = self.get_clean_fixed_surf()
            surf[:,:] = [[sit.what(res) for sit in ax] 
                for ax in self._lattice_[fixed]]
            return surf
        grids = []
        #for res in self.resources:
        for rsdx in range(len(self.resources)):
            grids.append(resource_grid(rsdx))
        return grids

def_string = ''
def simulate(sys_string = def_string):
    species = ['A', 'B', 'C']
    #agents = [
    #    ((1,3),'A'),((4,3),'B'),((8,7),'C'), 
    #    ((95,7),'A'),((84,9),'B'),((88,2),'C'), 
    #    ((93,83),'A'),((84,93),'B'),((88,87),'C'), 
    #    ((2,93),'A'),((6,83),'B'),((4,97),'C'), 
    #        ]
    agents =\
        [('random', 'A')]*3333 +\
        [('random', 'B')]*3333 +\
        [('random', 'C')]*3333 
    reactions = [
        #'(1)A+(1)food[5.0]->(2)A', 
        #'(1)B+(1)food[5.0]->(2)B', 
        #'(1)C+(1)food[5.0]->(2)C', 
        '(1)A+(1)B[5.0]->(2)A', 
        '(1)B+(1)C[5.0]->(2)B', 
        '(1)C+(1)A[5.0]->(2)C', 
        #'(1)A[0.01]->', 
        #'(1)B[0.01]->', 
        #'(1)C[0.01]->', 
        #'[0.01]->(1)A', 
        #'[0.01]->(1)B', 
        #'[0.01]->(1)C', 
                ]
    #birthing_rxns = [6]
    birthing_flag = False#True
    resources = []#[('all', 'food', 3)]
    lattice_dims = 2
    lattice_size = 100
    max_occupancy = 2
    timed_out = False
    timed_out_limit = 3600.0

    # I TREAT AGENTS AT EACH LATTICE SITE AS INDISTINGUISHABLE
    #  BUT THE LOCAL AGENTS ARE REMOVED IN REVERSE ORDER OF ENTERING THE SITE
    ### to properly include birth reactions - treat lattice sites as agents
    ###  with zero motion propensity
    # still want occupancy handling
    # species can only have a single letter to identify
    # voxel visualization pipeline for 3-d subspaces (perfectly analogous to
    #  surface_data pipeline)
    # statistics on average agent (average agent trajectory really....)
    # add full function/variable pipeline for reactions
    # add full mcfg/sys_string pipeline
    # introduce concept of nations (groups of agents)
    # i double count propensities since which agent initiates a reaction
    #  is never relevant - only the variety of agents post reaction

    time = 0.0
    last_time = 0.0
    end_time = 50.0
    incr_time = 0.1

    total_captures = end_time/incr_time
    capture_count = 0                               
    res_names = [r[1] for r in resources]
    targets = ['time', 'population', 'mean_occupancy'] + species[:]# + res_names
    plot_targets = ['time', 'population', 'mean_occupancy'] + species[:]# + res_names
    surf_targets = [
        'population_surfaces', 
        'identity_surfaces'] + ['resource_surface_' + r for r in res_names]
    plot_surf_targets = [
        'population_surfaces', 
        'identity_surfaces'] + ['resource_surface_' + r for r in res_names] 
    target_dexes = [targets.index(ta) for ta in plot_targets]
    surf_target_dexes = [surf_targets.index(ta) for ta in plot_surf_targets]
    _lattice_ = lattice(
        dims = lattice_dims, ax_length = lattice_size, species = species, 
        reactions = reactions, max_occupancy = max_occupancy, 
        birthing_flag = birthing_flag, resources = resources)
    for ag in agents: _lattice_.put_onto_lattice_at_site(ag[0],agent(ag[1]))
    for ag in _lattice_.agents: ag.update_propensity(_lattice_)
    maxax = _lattice_.max_axis_size
    state = np.zeros(shape = (len(targets)), dtype = np.float)
    surface_state = np.zeros(shape=(len(surf_targets),maxax,maxax),dtype=np.float)
    data = np.zeros(shape=(len(plot_targets),total_captures),dtype=np.float)
    surface_data = np.zeros(shape = (len(surf_targets), 
    	total_captures, maxax, maxax), dtype = np.float)
    set_surface_state(surface_state, time, _lattice_)
    set_state(state, surface_state, time, _lattice_)

    start_time = ti.time()
    if birthing_flag: dead_eco_check = lambda : False
    else: dead_eco_check = lambda : not _lattice_.agent_count > 0
    dead_eco = dead_eco_check()
    while capture_count < total_captures and not dead_eco and not timed_out:
        propensity_total = _lattice_.total_propensity()
        if propensity_total > 0.0 and _lattice_.agent_count > 0:
            propensity_total_inv = 1.0/propensity_total
            time_step = -1.0 * log(rm.random()) * propensity_total_inv
            acnt = len(_lattice_.agents)
            agrand = rm.uniform(0.0, propensity_total)
            agent_dex = pick_agent(_lattice_.agents, acnt, agrand)
        else:
            time_step = incr_time
            agent_dex = -1

        time += time_step
        #set_surface_state(surface_state, time, _lattice_)
        set_state(state, surface_state, time, _lattice_)

        real_time = state[0]
        if last_time < real_time and capture_count < total_captures:
            set_surface_state(surface_state, time, _lattice_)
            set_state(state, surface_state, time, _lattice_)
        while last_time < real_time and capture_count < total_captures:
            state[0] = last_time
            last_time += incr_time
            capture(data, state, capture_count, target_dexes)
            capture(surface_data, surface_state, capture_count, surf_target_dexes)
            capture_count += 1
            print 'time', last_time, 'pop', len(_lattice_.agents)
        state[0] = real_time

        #agent activity cant depend on state/surface_state...
        if agent_dex >= 0: _lattice_.act_agent(agent_dex)
        agents_propensity(_lattice_)

        if ti.time() - start_time > timed_out_limit: timed_out = True
        dead_eco = dead_eco_check()

    if timed_out or dead_eco:                        
        #set_surface_state(surface_state, time, _lattice_)
        state[0] = last_time
        last_time += incr_time
        capture(data, state, capture_count, target_dexes)
        capture(surface_data, surface_state, capture_count, surf_target_dexes)
        capture_count += 1
        print 'time', last_time, 'pop', len(_lattice_.agents)

        toss = capture_count
        surface_data = surface_data[:,:toss,:,:]
        data = data[:,:toss]
        return data,surface_data,targets + surf_targets
    return data, surface_data, targets + surf_targets

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1: use_pyx = sys.argv[1]
    else: use_pyx = False
    if use_pyx:
        import lattice_simulator as lsim
        print 'USING CYTHON!'
        lprf.profile_function(lsim.simulate)
    else:
        lprf.profile_function(simulate)
    #lprf.profile_function(simulate)
















