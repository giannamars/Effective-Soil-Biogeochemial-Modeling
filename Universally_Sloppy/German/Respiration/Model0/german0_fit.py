from SloppyCell.ReactionNetworks import *
import numpy as np
import scipy
import german_expt

model_net = IO.from_SBML_file('germans.xml', 'base')
model_net.set_var_optimizable('I', False)

# Define respiration
model_net.add_species('Resp', 'compartmentOne')
model_net.add_assignment_rule('Resp', '1 - eps*Vmax*x_1*x_2/(Km + x_1)')

## Output latex'ed equations for debugging
IO.eqns_TeX_file(model_net, 'model0.tex')


# Create the model
m = Model([german_expt.expt], [model_net])
p0= m.get_params().copy()
print(len(p0))


## Set prior ranges from value/prior_range to value*prior_range in interval [0,1] for substrate utilization efficiencies
res = Residuals.PriorInLog('eps_prior', 'eps', np.log(p0[0]), np.log(np.sqrt(p0[0])))
m.AddResidual(res)


## Optimize to fit data
print 'Initial Cost:', m.cost(p0)
popt = Optimization.fmin_lm_log_params(m, p0, maxiter=100, disp=True)
cost_opt = m.cost(popt)
print(popt)
print 'Optimized Cost:', m.cost(popt)
Utility.save(popt, 'popt0.model.bpkl')


## Ensemble fit
#print('Calculating Hessian')
#j = m.jacobian_log_params_sens(np.log(popt))
#jtj = np.dot(np.transpose(j),j)
#np.savetxt('hessian0.dat', jtj)
#e,v = Utility.eig(jtj)
#e = scipy.real(e)
#print('Generating ensemble')
#ens_data, costs, ratio = Ensembles.ensemble_log_params(m, popt, jtj, steps=1000000, recalc_hess_alg = True)
##
#ens_data = np.array(ens_data)
#Utility.save(ens_data, 'model0.ens_data_jac.bpkl')
#costs = np.array(costs)
#Utility.save(costs, 'model0.costs_jac.bpkl')


# We can look at the autocorrelation function of our costs to check efficiency
#Plotting.figure(1)
#ac = Ensembles.autocorrelation(costs)
#Plotting.plot(ac)
