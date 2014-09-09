#cython: profile=True,boundscheck=False,cdivision=True

import random as rm
import numpy as np





cpdef float propensity(list used, str populace, float rate):
    cdef int upop = 1
    cdef int tupop
    cdef int tucnt
    cdef int udex
    cdef int s
    cdef float prop
    for udex in xrange(len(used)):
        u = used[udex]
        if u[2]: continue
        tupop = 0
        tucnt = populace.count(u[1]) 
        for s in xrange(u[0]):
            tupop += tucnt - s
            upop *= tupop
    prop = upop * rate
    return prop

cpdef bint validate_rxn(list used, str pop, list resources):
    for u in used:
        if not u[2]:
            if not pop.count(u[1]) >= u[0]:
                return False
        else:
            rdex = [res.variety for res in resources].index(u[1])
            resource = resources[rdex]
            if not resource.quantity >= u[0]: return False
    return True

cpdef int choose_biased_index(list bins, int bcnt):
    cdef float tot = 0
    cdef int l
    cdef list lookup = []
    cdef float rand
    for b in xrange(bcnt):
        tot += bins[b]
        lookup.append(tot)

    if tot == 0.0:
        l = rm.randrange(bcnt)
        return l

    rand = rm.uniform(0.0, tot)
    for l in xrange(bcnt):
        if rand < lookup[l]:
            return l

#1000000 is a HARD limit on the maximum number of agents!!
cdef float prop_table[1000000]
def pick_agent_to_act(list agents, int pcnt, float rand):
    cdef float tot = 0.0
    cdef float agtot
    cdef int adx
    for adx in xrange(pcnt):
        agtot = agents[adx].total_prop
        tot += agtot
        prop_table[adx] = tot

    for adx in xrange(pcnt):
        if rand < prop_table[adx]:
            return adx

cpdef update_propensities(lattice):
    cdef int ldex
    cdef int adex
    cdef int fcnt
    cdef int scnt
    flagged = lattice.flagged_locs
    fcnt = len(flagged)
    latt = lattice._lattice_
    for ldex in xrange(fcnt):
        flag = flagged[ldex]
        sags = latt[flag].agents
        scnt = len(sags)
        for adex in xrange(scnt):
            sags[adex].update_propensity(lattice)
    lattice.flagged_locs = []

cpdef float update_agents_react_prop(site, list rxns, unsigned int rxcnt):
    #cdef int rxcnt = len(rxns)
    cdef float prop = 0.0
    populace = site.populace
    resources = site.resources
    for rdx in xrange(rxcnt):
        rx = rxns[rdx]
        used = rx.used
        if validate_rxn(used, populace, resources):
            rate = rx.rate
            prop += propensity(used, populace, rate)
    return prop






