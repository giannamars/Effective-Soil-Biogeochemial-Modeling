# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 12:22:00 2018

@author: Gianna
"""

import copy
from SloppyCell.ReactionNetworks import *
import numpy as np
import german_expt

model_net = IO.from_SBML_file('germans.xml', 'base')
model_net.set_var_optimizable('I', False)
# Define respiration
model_net.add_species('Resp', 'compartmentOne')
model_net.add_assignment_rule('Resp', '1 - eps*Vmax*x_1*x_2/(Km + x_1)')

# Create the model
m = Model([german_expt.expt], [model_net])
noprior_m = copy.deepcopy(m)

p0 = m.get_params().copy()
# This is a prior that runs (with 95% probability) from value/prior_range to


res = Residuals.PriorInLog('eps_prior', 'eps', np.log(p0[0]), np.log(np.sqrt(p0[0])))
m.AddResidual(res)

# Optimize to fit data
print 'Initial Cost:', m.cost(p0)
popt = Optimization.fmin_lm_log_params(m, p0, maxiter=100, disp=True)

import matplotlib.pyplot as plt
import matplotlib

fig_width = 3  # inches
golden_mean = (np.sqrt(5)-1.0)/2.0   # aethestic ratio
fig_height = fig_width*golden_mean
dpi = 300.0

fig_size = [fig_width, fig_height]

params = {'backend': 'ps',
          'axes.labelsize': 8,
          'axes.titlesize': 8,
          'font.size': 8,
          'legend.fontsize': 8,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'text.usetex': False,
          'figure.figsize': fig_size}


plt.rcParams.update(params)

fig, ax_traj = plt.subplots(figsize=fig_size, dpi=300)
fig.subplots_adjust(left=.18, bottom=.13, right=.99, top=.97)

traj = Dynamics.integrate(model_net, [0, 8640],rtol=1e-9, params=popt, fill_traj=True)
ax_traj.plot(traj.get_times(), traj.get_var_traj('x_1'), color='#525252', linewidth=2)
expt = m.exptColl.values()[0]
tt = expt.data['base']['x_1'].keys()
yy = [expt.data['base']['x_1'][t][0] for t in tt]
yerr = [expt.data['base']['x_1'][t][1] for t in tt]
ax_traj.errorbar(tt,yy,yerr,fmt='.', elinewidth=1,color='xkcd:black', ms=7, label='Perfect')

ax_traj.set_ylabel('S')
ax_traj.set_xlabel('Time')
fig.savefig('german_S_0_fit.pdf')

fig, ax_traj = plt.subplots(figsize=fig_size, dpi=300)
fig.subplots_adjust(left=.18, bottom=.13, right=.99, top=.97)

traj = Dynamics.integrate(model_net, [0, 8640],rtol=1e-9, params=popt, fill_traj=True)
ax_traj.plot(traj.get_times(), traj.get_var_traj('x_2'), color='#525252', linewidth=2)
expt = m.exptColl.values()[0]
tt = expt.data['base']['x_2'].keys()
yy = [expt.data['base']['x_2'][t][0] for t in tt]
yerr = [expt.data['base']['x_2'][t][1] for t in tt]
ax_traj.errorbar(tt,yy,yerr,fmt='.', elinewidth=1,color='xkcd:black', ms=7, label='Perfect')

ax_traj.set_ylabel('B')
ax_traj.set_xlabel('Time')
fig.savefig('german_B_0_fit.pdf')

fig, ax_traj = plt.subplots(figsize=fig_size, dpi=300)
fig.subplots_adjust(left=.18, bottom=.13, right=.99, top=.97)

traj = Dynamics.integrate(model_net, [0, 8640],rtol=1e-9, params=popt, fill_traj=True)
ax_traj.plot(traj.get_times(), traj.get_var_traj('Resp'), color='#525252', linewidth=2)
expt = m.exptColl.values()[0]
tt = expt.data['base']['Resp'].keys()
yy = [expt.data['base']['Resp'][t][0] for t in tt]
yerr = [expt.data['base']['Resp'][t][1] for t in tt]
ax_traj.errorbar(tt,yy,yerr,fmt='.', elinewidth=1,color='xkcd:black', ms=7, label='Perfect')

ax_traj.set_ylabel('Resp')
ax_traj.set_xlabel('Time')
fig.savefig('german_Resp_0_fit.pdf')
