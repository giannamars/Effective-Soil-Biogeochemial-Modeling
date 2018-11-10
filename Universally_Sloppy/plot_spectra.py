import os
import scipy
import scipy.io
import numpy as np

from SloppyCell.ReactionNetworks import *
import pylab

fig_width = 3.487  # inches
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


pylab.rcParams.update(params)



import Common

pylab.figure(1)
pylab.clf()
leftmargin = 0.175
bottommargin = 0.11
rightmargin = 0.01
topmargin= 0.03
pylab.axes([leftmargin,bottommargin, 1.0 - rightmargin-leftmargin, 1.0 - topmargin-bottommargin])

ax = pylab.gca()

for model_ii, (model,temp) in enumerate(Common.model_list):
    h = np.loadtxt(os.path.join(model, 'hessian.dat'))

    e,v = Utility.eig(h)
    e = scipy.real(e)

    width = (234/len(e))**0.25  * 0.25
    l = Plotting.plot_eigval_spectrum(e/max(e), offset=0.15+model_ii,ax=ax, widths=0.7, lw = 1.0)


e = np.loadtxt(os.path.join('PECCATPCA', 'hessian.dat'))
e = scipy.real(e)
l = Plotting.plot_eigval_spectrum(e/max(e), offset=0.15+5, ax=ax, widths=0.7, lw=1.0)


for ii in range(1, model_ii+2):
    ax.plot([ii, ii], [1e-24, 2], '-', lw=0.5, color='black')

ax.set_ylim(1e-24, 2)

# Add labels
ax.set_xticks(0.5 + scipy.arange(model_ii+2))
import string
xlabels = ['%s' % tup for tup in string.ascii_lowercase[:model_ii+1]]
ax.set_xticklabels(xlabels, fontsize=10, verticalalignment='bottom', rotation=0, horizontalalignment='center')

for l in ax.get_xticklabels():
    l.set_y(l.get_position()[1] - 0.04)

for l in ax.get_xticklines():
    l.set_visible(False)

ax.set_xlim(0, model_ii+2)

ax.set_ylabel('Normalized Hessian eigenvalues')

pylab.savefig('universal.svg', transparent=True)


pylab.show()
