#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 09:36:31 2018

@author: giannamarschmann
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Feb 01 08:30:17 2018

@author: Gianna
"""

import os

import scipy
import scipy.io
import numpy

from SloppyCell.ReactionNetworks import *
#import matplotlib.pyplot as plt 
#import matplotlib
#matplotlib.rc('font',**{'family':'sans-serif',
#                        'sans-serif':['Arial'],
#                        'style':'normal',
#                        'size':12 })


import pylab

fig_width = 3.487  # inches
golden_mean = (numpy.sqrt(5)-1.0)/2.0   # aethestic ratio
fig_height = fig_width*golden_mean
#fig_height = 10.0
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

hblue = (0, 0.247, 0.459)
tred = (0.647, 0.118, 0.216)
#Plotting.figure(figsize=(6, 3.71))


f, (a1, a0) = pylab.subplots(2,1, gridspec_kw = {'height_ratios':[1, 5]})

#pylab.figure(1)
#pylab.clf()
#leftmargin = 0.15
#bottommargin = 0.1
#rightmargin = 0.07
#topmargin= 0.1
#pylab.axes([leftmargin,bottommargin, 1.0 - rightmargin-leftmargin, 1.0 - topmargin-bottommargin])
f.subplots_adjust(left=.18, bottom=.165, right=.99, top=.98)

for model_ii, (model,temp) in enumerate(Common.model_list):
    # Load the hessian
    model_iii = model_ii + 40
    h = numpy.loadtxt(os.path.join(model, 'hessian%s.dat' % model_iii))

    e, v = Utility.eig(h)
    e = scipy.real(e)

    width = (234/len(e))**0.25 * 0.25
    l = Plotting.plot_eigval_spectrum(e, offset = 0.15+model_ii, ax=a0, 
                                      widths=0.7, lw=1.0, lc='#cccccc')
    
# Now a lot of fiddling to make the plot prettier
#ax = Plotting.gca()
#for ii in range(1, model_ii+1):
#    ax.plot([ii, ii], [0.5e-27, 1], '-', lw=1, color="black")
a0.plot([0, 38], [1, 1], '--', lw=1, color='#ff6600')
    
# Add labels
a0.set_xticks(0.5 + scipy.arange(model_ii+1))
import string
xlabels = ['%s' % n for n in  numpy.append(scipy.arange(21,17,-1), [37,36])]
#xlabels=[]
a0.set_xticklabels(xlabels, verticalalignment='bottom',
                   rotation=0, horizontalalignment='center')
for l in a0.get_xticklines():
    l.set_visible(False)
a0.set_xlim(0, model_ii+1)
for l in a0.get_xticklabels():
    l.set_y(l.get_position()[1] - 0.055)


a0.set_ylim(1e-20, 1e9)
a0.set_ylabel('Hessian eigenvalues')
a0.set_xlabel('Number of parameters')


Jvals = numpy.array([1.93, 2.70, 3.11, 3.11, 15.12])
nums = numpy.arange(23)
pnums = nums[22:17:-1] 

a1.plot(pnums, Jvals, lw=1, color = 'k')
a1.set_xlim(21,18)
a1.set_xticklabels([])
a1.set_xticks([])

a1.set_ylim(1.5, 17)
a1.set_ylabel(r'$J(p)$')
a1.set_yticks([6,12])
#ax.text(2, 1e-10, 'Fungal co-metabolism', color=tred)
#ax.text(8, 1e-7, 'Maintenance respiration', color=tred)
#ax.text(14, 1e-4, 'Exponential microbial death', color=tred)

pylab.subplots_adjust(hspace=0.1)

f.savefig('reductionspectraplot2.svg',dpi=dpi, transparent=True)
pylab.show()
