# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 11:48:15 2018

@author: Gianna
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy
from SloppyCell.ReactionNetworks import *
from geodesic import geodesic, InitialVelocity

import model35_fit

x = np.log(model35_fit.popt)

# Calculate jtj approximation to hessian (in log params) and plot eigenvalue
# spectrum
j = model35_fit.m.jacobian_log_params_sens(np.log(model35_fit.popt))
jtj = np.dot(np.transpose(j), j)
np.savetxt('hessian35.dat', jtj)
e, v = Utility.eig(jtj)
e = scipy.real(e)
Plotting.figure(1)
l = Plotting.plot_eigval_spectrum(e, offset = 0.15, widths=0.7, lw=2)

## Reduce

func = lambda logp: np.array(model35_fit.m.res_log_params(logp))
jacobian = lambda logp: np.array(model35_fit.m.jacobian_log_params_sens(logp).values())

M = jacobian(x).shape[0]
N = jacobian(x).shape[1]
print(M)
print(N)
# Directional 2nd derivative
def Avv(x,v):
    h = 1e-2
    return (func(x + h*v) + func(x - h*v) - 2*func(x))/h/h

# Callback function used to monitor the geodesic after each step, integrate
# until norm of parameter velocity has grown by factor of 25
def callback(geo):
    # Integrate until the norm of the velocity has grown by a factor of 10
    # and print out some diagnotistic along the way
    print("Iteration: %i, tau: %f, |v| = %f" %(len(geo.vs), geo.ts[-1], np.linalg.norm(geo.vs[-1])))
    return np.linalg.norm(geo.vs[-1]) < 22.0*np.linalg.norm(geo.vs[0])


# Calculate initial velocity
v = InitialVelocity(x, jacobian, Avv)

# Calculate geodesic
geo = geodesic(func, jacobian, Avv, M, N, x, v,
                            atol = 1e-2, rtol = 1e-2, callback=callback, invSVD=True)
geo.integrate(15.0)

print('Got to t={0}'.format(geo.ts[-1]))
Plotting.figure(2)
plt.plot(geo.ts, geo.xs)

# Calculate jtj approximation to hessian (in log params) and plot eigenvalue
# spectrum
j = model35_fit.m.jacobian_log_params_sens(geo.xs[-1,:])
jtj = np.dot(np.transpose(j), j)
e, v = Utility.eig(jtj)
e = scipy.real(e)
Plotting.figure(3)
l = Plotting.plot_eigval_spectrum(e, offset = 0.15, widths=0.7, lw=2)


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

plt.show()
print('corresponding to %d',model35_fit.popt.items()[ SDisc.items()[0][0]])

# Save parameter values at end of geodesic
p0 = np.exp(geo.xs[-1,:])
Utility.save(p0, 'p35.geodesics.bpkl')

