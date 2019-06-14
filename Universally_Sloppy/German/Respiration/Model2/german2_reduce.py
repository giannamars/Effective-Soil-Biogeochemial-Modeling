# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 11:48:15 2018

@author: Gianna
"""

import numpy as np
import scipy
from SloppyCell.ReactionNetworks import *
from geodesic import geodesic, InitialVelocity

import matplotlib.pyplot as plt
import matplotlib

fig_width = 5  # inches
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

import german2_fit

x = np.log(german2_fit.popt)

# Calculate jtj approximation to hessian (in log params) and plot eigenvalue
# spectrum at the start of the geodesic
j = german2_fit.m.jacobian_log_params_sens(np.log(german2_fit.popt))
jtj = np.dot(np.transpose(j), j)
np.savetxt('hessian2.dat', jtj)
e, v = Utility.eig(jtj)
e = scipy.real(e)
Plotting.figure(1)
l = Plotting.plot_eigval_spectrum(e, offset = 0.15, widths=0.7, lw=2)


## Define residual and jacobian as function of log parameters
func = lambda logp: np.array(german2_fit.m.res_log_params(logp))
jacobian = lambda logp: np.array(german2_fit.m.jacobian_log_params_sens(logp).values())

M = jacobian(x).shape[0]
N = jacobian(x).shape[1]


## Define directional 2nd derivative
def Avv(x,v):
    h = 1e-1
    return (func(x + h*v) + func(x - h*v) - 2*func(x))/h/h


## Callback function used to monitor the geodesic after each step, integrate
# until norm of parameter velocity has grown by factor of 25
def callback(geo):
    # Integrate until the norm of the velocity has grown by a factor of 25
    # and print out some diagnotistic along the way
    print("Iteration: %i, tau: %f, |v| = %f" %(len(geo.vs), geo.ts[-1], np.linalg.norm(geo.vs[-1])))
    return np.linalg.norm(geo.vs[-1]) < 25.0*np.linalg.norm(geo.vs[0])


## Calculate initial velocity
v = InitialVelocity(x, jacobian, Avv)

## Calculate geodesic
geo = geodesic(func, jacobian, Avv, M, N, x, v,
                            atol = 1e-2, rtol = 1e-2, callback=callback)
geo.integrate(5.0)

print('Got to t={0}'.format(geo.ts[-1]))


fig, ax = plt.subplots(figsize=fig_size, dpi=300)
fig.subplots_adjust(left=.13, bottom=.13, right=.99, top=.85)

labels = ['Y', r'$\vartheta_1$']
colors = ['#cccccc', '#969696']
for i in range(geo.xs.shape[1]):
    ax.plot(geo.ts, geo.xs[:,i], label=labels[i],color=colors[i])
    
ax.set_xlabel(r'$\tau$')
ax.set_ylabel(r'$\log\, p$')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.16),
          ncol=4)

fig.savefig('german_limit_2.png', dpi=300)

## Calculate jtj approximation to hessian (in log params) and plot eigenvalue
# spectrum at end of geodesic
j = german2_fit.m.jacobian_log_params_sens(geo.xs[-1,:])
jtj = np.dot(np.transpose(j), j)
e, v = Utility.eig(jtj)
e = scipy.real(e)
Plotting.figure(3)
l = Plotting.plot_eigval_spectrum(e, offset = 0.15, widths=0.7, lw=2)


## Discern whether limit has been reached
def SDiscern(vec):
    Disc={}
    vnorm = vec/np.linalg.norm(vec)

    n=-1
    arglist = list(np.argsort(np.abs(vnorm)))
    arglist.reverse()
    for i in range(len(vnorm)-1):
        if np.abs(vnorm[arglist[i]]) - np.abs(vnorm[arglist[i+1]]) > 0.8/np.sqrt(i+1):
            n = i+1
            break
    if n == -1:
        return Disc

    for j in arglist[:n]:
        if vnorm[j]>0:
            Disc[j]=np.inf
        else:
            Disc[j] = 0

    return Disc



SDisc = SDiscern(geo.vs[-1]-geo.vs[-2])

print(SDisc)
print('corresponding to %d',german2_fit.popt.items()[ SDisc.items()[0][0]])

# Save parameter values at end of geodesic
p0 = np.exp(geo.xs[-1,:])
Utility.save(p0, 'p2.geodesics.bpkl')

plt.show()
