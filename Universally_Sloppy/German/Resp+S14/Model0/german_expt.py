# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 12:19:16 2018

@author: Gianna
"""

import numpy
from SloppyCell.ReactionNetworks import *

model_net = IO.from_SBML_file('germans_isotopes.xml', 'base')
model_net.set_var_optimizable('lambd', False)

# Define respiration
model_net.add_species('Resp', 'compartmentOne')
model_net.add_assignment_rule('Resp', '1 - eps*Vmax*x_1*x_2/(Km + x_1)')
# Define Delta14C
model_net.add_species('SOD14', 'compartmentOne')
model_net.add_assignment_rule('SOD14', '(x_3/x_1 - 1)*1000')
model_net.add_species('BOD14', 'compartmentOne')
model_net.add_assignment_rule('BOD14', '(x_4/x_2 - 1)*1000')

traj_full = Dynamics.integrate(model_net, [0,8640], fill_traj=True)

sample_traj_resp = Dynamics.integrate(model_net, numpy.linspace(0,8640,500), 
                                 fill_traj=False)

sample_traj_r14 = Dynamics.integrate(model_net, numpy.linspace(0,8640,12), 
                                 fill_traj=False)

sigma = 0.001

numpy.random.seed(290137)


times_resp = sample_traj_resp.get_times()[1:]
times_r14 = sample_traj_r14.get_times()[1:]
times_s14 = numpy.array([traj_full.get_times()[1],traj_full.get_times()[-1]])

data = {}
for var in ['Resp']:
    maxy = traj_full.get_var_traj(var).max()

    yy = sample_traj_resp.get_var_traj(var)[1:]
    noisy_yy = yy*(1+sigma*numpy.random.randn(len(times_resp)))
    uncert = sigma * maxy

    var_data = dict((t, (y, uncert)) for (t,y) in zip(times_resp, noisy_yy))
    data[var] = var_data
    
    
#for var in ['BOD14']:
#    maxy = traj_full.get_var_traj(var).max()
#
#    yy = sample_traj_r14.get_var_traj(var)[1:]
#    noisy_yy = yy*(1+sigma*numpy.random.randn(len(times_r14)))
#    uncert = sigma * maxy
#
#    var_data = dict((t, (y, uncert)) for (t,y) in zip(times_r14, noisy_yy))
#    data[var] = var_data
#    
    
for var in ['SOD14']:
    maxy = traj_full.get_var_traj(var).max()

    yy = numpy.array([traj_full.get_var_traj(var)[1], traj_full.get_var_traj(var)[-1]])
    noisy_yy = yy*(1+sigma*numpy.random.randn(len(times_s14)))
    uncert = sigma * maxy

    var_data = dict((t, (y, uncert)) for (t,y) in zip(times_s14, noisy_yy))
    data[var] = var_data
    
    
#for var in ['SOD14']:
#    maxy = traj_full.get_var_traj(var).max()
#
#    yy = traj_full.get_var_traj(var)[1:]
#    noisy_yy = yy*(1+sigma*numpy.random.randn(1))
#    uncert = sigma * maxy
#
#    var_data = {times_s14: (noisy_yy, uncert)}
#    data[var] = var_data
    
    

expt = Experiment('german_expt',
                  {'base': data}, 
                  fixedScaleFactors={'x_1':1,'x_2':1, 'Resp':1, 'BOD14':1, 'SOD14':1})