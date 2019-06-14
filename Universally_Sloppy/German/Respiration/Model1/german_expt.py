# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 12:19:16 2018

@author: Gianna
"""

import numpy
from SloppyCell.ReactionNetworks import *

model_net = IO.from_SBML_file('germans.xml', 'base')
model_net.set_var_optimizable('I', False)

# Define respiration
model_net.add_species('Resp', 'compartmentOne')
model_net.add_assignment_rule('Resp', '1 - eps*Vmax*x_1*x_2/(Km + x_1)')

traj_full = Dynamics.integrate(model_net, [0,8640], fill_traj=True)

sample_traj = Dynamics.integrate(model_net, numpy.linspace(0,8640,500), 
                                 fill_traj=False)

sigma = 0.001

numpy.random.seed(290137)


times = sample_traj.get_times()[1:]
data = {}
for var in ['Resp']:
    maxy = traj_full.get_var_traj(var).max()

    yy = sample_traj.get_var_traj(var)[1:]
    noisy_yy = yy*(1+sigma*numpy.random.randn(len(times)))
    uncert = sigma * maxy

    var_data = dict((t, (y, uncert)) for (t,y) in zip(times, noisy_yy))
    data[var] = var_data
    
    

expt = Experiment('german_expt',
                  {'base': data}, 
                  fixedScaleFactors={'x_1':1,'x_2':1, 'Resp':1})