# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 12:55:17 2018

@author: Gianna
"""

from numpy import *
import scipy.stats
from SloppyCell.ReactionNetworks import *
import monod_net, monod_model
import matplotlib.pyplot as plt

# Calculate an ensemble
#

print("Generating ensemble")
Network.full_speed()
ens_data, costs, ratio = Ensembles.ensemble_log_params(monod_model.m,
                                                       monod_model.popt,
                                                       hess=monod_model.jtj,
                                                       steps=100,
                                                       seeds=(113,207),
                                                       recalc_hess_alg=False)

ens_data = array(ens_data)
Utility.save(ens_data, 'monod.ens_data.bpkl')

costs = array(costs)
Utility.save(costs, 'monod.costs.bpkl')

print "Ensemble has %i members. Acceptance ratio was %f." % (len(ens_data), ratio)


# Calcualate cost surface
print("Generating cost surface")
Npt = 151
land_xx = logspace(-3, 2, Npt)
land_yy = logspace(-3, 2, Npt)
Z = zeros((len(land_yy), len(land_xx)))
for ii,y in enumerate(land_yy):
    for jj,x in enumerate(land_xx):
        p = [x,y]
        Z[ii,jj] = monod_model.noprior_m.cost(p)
Utility.save((land_xx, land_yy, Z), 'monod.model_surface.bpkl')
