from SloppyCell.ReactionNetworks import *
import PECCAT_experiment_MCPA
import numpy as np
import scipy

## Read model from SBML file
model_net = IO.from_SBML_file('PECCAT32.xml', 'base')

## define bulk measurements
model_net.add_species('cmic', 'compartmentOne')
model_net.add_assignment_rule('cmic', 'x_1 + x_2 + x_3')
#
model_net.add_species('cBtot', 'compartmentOne')
model_net.add_assignment_rule('cBtot', 'x_1 + x_2')
#
model_net.add_species('doc', 'compartmentOne')
model_net.add_assignment_rule('doc', 'x_6 + x_8')
#
model_net.add_species('toc', 'compartmentOne')
model_net.add_assignment_rule('toc', 'x_1 + x_2 + x_3 + x_6 + x_7 + x_8 + x_9')

## Add inital activity fractions to parameters
model_net.add_parameter('rF0', 0.05)
#
model_net.set_var_ic('x_5', 'rF0')
#
model_net.add_assignment_rule('x_4', ' Phi21*x_8')
#

## Don't optimize MCPA sorption kinetics and soil properties
model_net.set_var_optimizable('Kf', False)
model_net.set_var_optimizable('nf', False)
model_net.set_var_optimizable('nc', False)
model_net.set_var_optimizable('theta', False)
model_net.set_var_optimizable('rhoB', False)
model_net.set_var_optimizable('cL', False)
model_net.set_var_optimizable('bL', False)
model_net.set_var_optimizable('YLHiq', True)
model_net.set_var_optimizable('YLLoq', True)


## Output latex'ed equations for debugging
IO.eqns_TeX_file(model_net, 'model32.tex')

## Create the model
m = Model([PECCAT_experiment_MCPA.expt], [model_net])
p0= m.get_params().copy()

#
p0 = Utility.load('poptnew1.model.bpkl')
#
p0[5]=100
p0[22]=0.95


## Set prior ranges from value/prior_range to value*prior_range

prior_range = 100
for key, value in p0.items()[0:len(p0)-10]:
    res = Residuals.PriorInLog('%s_prior' % key, key, np.log(value),
                               np.log(np.sqrt(prior_range)))
    m.AddResidual(res)


res = Residuals.PriorInLog('YLHiq_prior', 'YLHiq', np.log(p0[14]), np.log(np.sqrt(p0[14])))
m.AddResidual(res)
res = Residuals.PriorInLog('YLLoq_prior', 'YLLoq', np.log(p0[15]), np.log(np.sqrt(p0[15])))
m.AddResidual(res)
res = Residuals.PriorInLog('YrB_prior', 'YrB', np.log(p0[16]), np.log(np.sqrt(p0[16])))
m.AddResidual(res)
res = Residuals.PriorInLog('YrF_prior', 'YrF', np.log(p0[17]), np.log(np.sqrt(p0[17])))
m.AddResidual(res)
res = Residuals.PriorInLog('YRFP_prior', 'YRFP', np.log(p0[18]), np.log(np.sqrt(p0[18])))
m.AddResidual(res)
res = Residuals.PriorInLog('YsBHiq_prior', 'YsBHiq', np.log(p0[19]), np.log(np.sqrt(p0[19])))
m.AddResidual(res)#
res = Residuals.PriorInLog('YsBLoq_prior', 'YsBLoq', np.log(p0[20]), np.log(np.sqrt(p0[20])))
m.AddResidual(res)
res = Residuals.PriorInLog('YsBPLoq_prior', 'YsBPLoq', np.log(p0[21]), np.log(np.sqrt(p0[21])))
m.AddResidual(res)
res = Residuals.PriorInLog('YsFLoq_prior', 'YsFLoq', np.log(p0[22]), np.log(np.sqrt(p0[22])))
m.AddResidual(res)
res = Residuals.PriorInLog('rF0_prior', 'rF0', np.log(p0[23]), np.log(np.sqrt(p0[23])))
m.AddResidual(res)

## Optimize to fit data
print 'Initial Cost:', m.cost(p0)
popt = Optimization.fmin_lm_log_params(m, p0, maxiter=1, disp=True)
# Then we run Levenberg-Marquardt
#popt = Optimization.leastsq_log_params(m, popt1)
cost_opt = m.cost(popt)
print(popt)
print(len(popt))
print 'Optimized Cost:', m.cost(popt)

Utility.save(popt, 'poptnew2.model.bpkl')



## Ensemble fit
print('Calculating Hessian')
j = m.jacobian_log_params_sens(np.log(popt))
jtj = np.dot(np.transpose(j),j)

#print('Generating ensemble')
#Network.full_speed()
#ens_data, costs, ratio = Ensembles.ensemble_log_params(m, popt, jtj, steps=100000)
#
#ens_data = np.array(ens_data)
#Utility.save(ens_data, 'model320.ens_data1.bpkl')
#costs = np.array(costs)
#Utility.save(costs, 'model320.costs.bpkl')

# We can look at the autocorrelation function of our costs to check efficiency
#Plotting.figure(1)
#ac = Ensembles.autocorrelation(costs)
#Plotting.plot(ac)


