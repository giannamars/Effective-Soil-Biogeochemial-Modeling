# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 12:22:00 2018

@author: Gianna
"""

import copy
from SloppyCell.ReactionNetworks import *
import numpy
import monod_expt, monod_net

# Create the model
m = Model([monod_expt.expt], [monod_net.net])
noprior_m = copy.deepcopy(m)

p0 = m.get_params().copy()
# This is a prior that runs (with 95% probability) from value/prior_range to
#  value*prior_range
prior_range = 1e3
for key, value in p0.items():
    res = Residuals.PriorInLog('%s_prior' % key, key, numpy.log(value),
                               numpy.log(numpy.sqrt(prior_range)))
    m.AddResidual(res)

# Optimize to fit the synthetic data set
popt = Optimization.leastsq_log_params(m, p0)
cost_opt = m.cost(popt)
print('Optimized cost is {}'.format(cost_opt))

print(popt)

# Calculate jtj approximation to hessian (in log params)
j = m.jacobian_log_params_sens(numpy.log(popt))
jtj = numpy.dot(numpy.transpose(j), j)
numpy.savetxt('monod.jtj.popt.dat', jtj)

e,v = Utility.eig(jtj)
print(e)

# Calculate jtj approximation to hessian (in log params)
j_noprior = noprior_m.jacobian_log_params_sens(numpy.log(popt))
jtj_noprior = numpy.dot(numpy.transpose(j_noprior), j_noprior)
numpy.savetxt('monod.jtj.popt.noprior.dat', jtj_noprior)

jtj_uncerts = numpy.sqrt(numpy.diag(numpy.linalg.inv(jtj)))

