# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 12:53:14 2018

@author: Gianna
"""

"""
Plots Fig. 2
"""
from SloppyCell.ReactionNetworks import *
import scipy.stats
from numpy import *
import monod_model
import matplotlib

import matplotlib.pyplot as plt


plt.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
plt.rc('text', usetex=False)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', labelsize=8)
plt.rc('legend', fontsize=8)
width = 3.0
golden_mean = (sqrt(5)-1.0)/2.0   # aethestic ratio
height = width*golden_mean


fig, ax_land = plt.subplots()
fig.subplots_adjust(left=.16, bottom=.225, right=.995, top=.96)



# Landscape plot
land_xx, land_yy, Z = Utility.load('monod.model_surface.bpkl')
Z[isnan(Z)] = Z[logical_not(isnan(Z))].max()


norm = matplotlib.colors.LogNorm(vmin=monod_model.cost_opt, vmax=Z.max())
ens_data = Utility.load('monod.ens_data.bpkl')
mappable = ax_land.pcolor(land_xx, land_yy, Z, norm=norm, cmap='Greys',antialiased=False)
ax_land.set_xscale('log')
ax_land.set_yscale('log')
ax_land.set_xlim(1e-3, 1e2)
ax_land.set_ylim(1e-3, 1e2)


# Colorbar
cbar = fig.colorbar(mappable)
cbar.set_ticks(concatenate([[monod_model.cost_opt],[10,100,1000]]))
cbar.set_ticklabels([r'${0:.1f}$'.format(monod_model.cost_opt)] + ['10', '100', '1000'] )

## Plot 100 ensemble fits
ens_data = Utility.load('monod.ens_data.bpkl')

def prune_ens(ens, Npt):
    return ens[1::len(ens)//Npt]


ens_data_pruned = prune_ens(ens_data, 120)
ax_land.plot(ens_data_pruned[:,0], ens_data_pruned[:,1], 'o', mfc='w', mec='k', ms=3, mew=0.1,  zorder=1)


#Plot  Best fit
ax_land.scatter([monod_model.popt[0]],[ monod_model.popt[1]], c='#fb8072', s=16, edgecolor='k', lw=0.1, zorder=50)


# Plot geodesics
geodesics = Utility.load('monod.geodesics.bpkl')
colors = iter(['#bebada', '#80b1d3'])
for xs, vs, ts in geodesics:
    ax_land.plot(exp(xs[:,0]), exp(xs[:,1]), c=next(colors), lw=1, zorder=2)

ax_land.set_xlabel(r'$\log\, V_{max}$')
ax_land.set_ylabel(r'$\log\, K_{S}$')
ax_land.xaxis.labelpad = 2
ax_land.yaxis.labelpad = 1


# Plot Hessian approximation

def contour_top_bottom(h, delC=1, Npts=100):
    """
    Contours tracing ellipse from the hessian approximation
    Note that these are countours in delta X and delta Y.
    """

    delx = sqrt(delC*h[1,1]/(h[1,1]*h[0,0] - h[0,1]**2))
    xx = linspace(-delx, delx, Npts)
    a = h[1,1]
    b = 2*h[0,1]*xx
    c = h[0,0]*xx**2 - delC
    # Avoid errors due to tiny negative values
    inner = b**2 - 4*a*c
    inner = maximum(inner, 0)
    top = (-b + sqrt(inner))/(2*a)
    bottom = (-b - sqrt(inner))/(2*a)

    return xx, top, bottom


# Hessian 95% confidence interval
delxx,delyy_top,delyy_bottom = contour_top_bottom(monod_model.jtj,
                                                  delC=4, Npts=100)
# Convert from deltas in log params to normal params
xx = exp(log(monod_model.popt[0]) + delxx)
yy_top = exp(log(monod_model.popt[1]) + delyy_top)
yy_bottom = exp(log(monod_model.popt[1]) + delyy_bottom)
# This concatenation ensures we plot a single smooth curve.
xx = concatenate((xx, xx[::-1], [xx[0]]))
yy = concatenate((yy_top, yy_bottom[::-1], [yy_top[0]]))
ax_land.plot(xx, yy, '-', c='#8dd3c7', lw=1, zorder=3)


# Plot other three points
ax_land.plot(2.5, 0.2, 'o', mfc='#fdb462', mec='k', ms=4, mew=0.1,  zorder=1)
ax_land.plot(5, 1.1, 'o', mfc='#b3de69', mec='k', ms=4, mew=0.1,  zorder=1)
ax_land.plot(0.15, 0.3, 'o', mfc='#ffffb3', mec='k', ms=4, mew=0.1,  zorder=1)



plt.minorticks_off()
ax_land.grid(False)


fig.set_size_inches(width, height)
fig.savefig('cost_landscape.png', dpi=300)

plt.show()
