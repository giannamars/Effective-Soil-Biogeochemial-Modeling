from SloppyCell.ReactionNetworks import *
import PECCAT_experiment
import numpy as np
import scipy

## Read model from SBML file
model_net = IO.from_SBML_file('PECCAT3.xml', 'base')

## define bulk measurements
model_net.add_species('cmic', 'compartmentOne')
model_net.add_assignment_rule('cmic', 'x_1 + x_2 + x_3')
#
model_net.add_species('cBtot', 'compartmentOne')
model_net.add_assignment_rule('cBtot', 'x_1 + x_2')
#
model_net.add_species('doc', 'compartmentOne')
model_net.add_assignment_rule('doc', 'x_7 + x_9')
#
model_net.add_species('toc', 'compartmentOne')
model_net.add_assignment_rule('toc', 'x_1 + x_2 + x_3 + x_7 + x_8 + x_9 + x_10')

## Add inital activity fractions to parameters
model_net.add_parameter('rB0', 0.1)
model_net.add_parameter('rBP0', 0.2)
model_net.add_parameter('rF0', 0.1)
#
model_net.set_var_ic('x_4', 'rB0')
model_net.set_var_ic('x_5', 'rBP0')
model_net.set_var_ic('x_6', 'rF0')
#
model_net.add_species('co2_tot', 'compartmentOne')
model_net.add_assignment_rule('co2_tot', 'x_11 + x_12*(1-YLHiq*(((cL + (t/(t**2 + bL))**3) - cL)/(cL + (1/2*sqrt(bL))**3)) - YLLoq*(1-(((cL + (t/(t**2 + bL))**3) - cL)/(cL + (1/2*sqrt(bL))**3))))*(cL + (t/(t**2 + bL))**3)')

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
IO.eqns_TeX_file(model_net, 'model3.tex')

## Create the model
m = Model([PECCAT_experiment.expt], [model_net])
p0= m.get_params().copy()

p3 = Utility.load('p2.geodesics.bpkl')
p3 = np.delete(p3, 29)


for i, (key,value) in enumerate(p0.items()):
    p0[i] = p3[i]

## Set prior ranges from value/prior_range to value*prior_range


res = Residuals.PriorInLog('YLHiq_prior', 'YLHiq', np.log(p0[41]), np.log(np.sqrt(p0[41])))
m.AddResidual(res)
res = Residuals.PriorInLog('YLLoq_prior', 'YLLoq', np.log(p0[42]), np.log(np.sqrt(p0[42])))
m.AddResidual(res)
res = Residuals.PriorInLog('YrB_prior', 'YrB', np.log(p0[43]), np.log(np.sqrt(p0[43])))
m.AddResidual(res)
res = Residuals.PriorInLog('YrF_prior', 'YrF', np.log(p0[44]), np.log(np.sqrt(p0[44])))
m.AddResidual(res)
res = Residuals.PriorInLog('YRFP_prior', 'YRFP', np.log(p0[45]), np.log(np.sqrt(p0[45])))
m.AddResidual(res)
res = Residuals.PriorInLog('YsBHiq_prior', 'YsBHiq', np.log(p0[46]), np.log(np.sqrt(p0[46])))
m.AddResidual(res)#
res = Residuals.PriorInLog('YsBLoq_prior', 'YsBLoq', np.log(p0[47]), np.log(np.sqrt(p0[47])))
m.AddResidual(res)
res = Residuals.PriorInLog('YsBPHiq_prior', 'YsBPHiq', np.log(p0[48]), np.log(np.sqrt(p0[48])))
m.AddResidual(res)
res = Residuals.PriorInLog('YsBPLoq_prior', 'YsBPLoq', np.log(p0[49]), np.log(np.sqrt(p0[49])))
m.AddResidual(res)
res = Residuals.PriorInLog('YsBPP_prior', 'YsBPP', np.log(p0[50]), np.log(np.sqrt(p0[50])))
m.AddResidual(res)
res = Residuals.PriorInLog('YsFHiq_prior', 'YsFHiq', np.log(p0[51]), np.log(np.sqrt(p0[51])))
m.AddResidual(res)
res = Residuals.PriorInLog('YsFLoq_prior', 'YsFLoq', np.log(p0[52]), np.log(np.sqrt(p0[52])))
m.AddResidual(res)
res = Residuals.PriorInLog('rB0_prior', 'rB0', np.log(p0[53]), np.log(np.sqrt(p0[53])))
m.AddResidual(res)
res = Residuals.PriorInLog('rBP0_prior', 'rBP0', np.log(p0[54]), np.log(np.sqrt(p0[54])))
m.AddResidual(res)
res = Residuals.PriorInLog('rF0_prior', 'rF0', np.log(p0[55]), np.log(np.sqrt(p0[55])))
m.AddResidual(res)

## Optimize to fit data
print 'Initial Cost:', m.cost(p0)
popt = Optimization.fmin_lm_log_params(m, p0, maxiter=5, disp=True)
# Then we run Levenberg-Marquardt
#popt = Optimization.leastsq_log_params(m, popt1)
cost_opt = m.cost(popt)
print(popt)
print(len(popt))
print 'Optimized Cost:', m.cost(popt)

Utility.save(popt, 'popt3.model.bpkl')
