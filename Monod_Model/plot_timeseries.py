# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 17:56:27 2018

@author: Gianna
"""

"""
Plots Fig. 2
"""
from SloppyCell.ReactionNetworks import *
import scipy.stats
from numpy import *
import monod_model

import matplotlib.pyplot as plt

plt.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
plt.rc('text', usetex=False)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', labelsize=8)
width = 3.0
golden_mean = (sqrt(5)-1.0)/2.0   # aethestic ratio
height = width*golden_mean



fig, ax_traj = plt.subplots()
fig.subplots_adjust(left=.18, bottom=.19, right=.96, top=.97)


net = monod_model.m.calcColl[0]

Npt = 31
land_xx = logspace(-3, 2, Npt)
land_yy = logspace(-3, 2, Npt)

# Plot family of trajectories

for ii,y in enumerate(land_yy):
    for jj,x in enumerate(land_xx):
        p = [x, y]
        traj = Dynamics.integrate(net, [0, 40], params=p)
        ax_traj.plot(traj.get_times(), traj.get_var_traj('x_2'), 'lightgray', linewidth=0.05)


# Highlight three special trajectories


p = [2.5, 0.2]
traj = Dynamics.integrate(net, [0, 40], params=p)
ax_traj.plot(traj.get_times(), traj.get_var_traj('x_2'), '#fdb462', linewidth=1)

p = [5,1.1]
traj = Dynamics.integrate(net, [0, 40], params=p)
ax_traj.plot(traj.get_times(), traj.get_var_traj('x_2'), '#b3de69', linewidth=1)

p = [0.15, 0.3]
traj = Dynamics.integrate(net, [0, 40], params=p)
ax_traj.plot(traj.get_times(), traj.get_var_traj('x_2'), '#ffffb3', linewidth=1)


traj = Dynamics.integrate(net, [0, 40], params=monod_model.popt)
ax_traj.plot(traj.get_times(), traj.get_var_traj('x_2'), '#fb8072', linewidth=1)
expt = monod_model.m.exptColl.values()[0]
tt = expt.data['monod']['x_2'].keys()
yy = [expt.data['monod']['x_2'][t][0] for t in tt]
yerr = [expt.data['monod']['x_2'][t][1] for t in tt]
ax_traj.errorbar(tt, yy, yerr, fmt='*',markersize=3, elinewidth=1, color='xkcd:black')
ax_traj.set_xlim(0, 40)
ax_traj.set_ylim(0.03, 0.285)

ax_traj.set_xlabel('Time [h]')
ax_traj.set_ylabel('Biomass [mg C/g]')


fig.set_size_inches(width, height)
fig.savefig('timeseries.pdf')


plt.show()




