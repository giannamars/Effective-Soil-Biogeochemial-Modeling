from SloppyCell.ReactionNetworks import *
import PECCAT_experiment
import numpy as np
import scipy

## Read model from SBML file
model_net = IO.from_SBML_file('PECCAT36.xml', 'base')

## define bulk measurements
model_net.add_species('cmic', 'compartmentOne')
model_net.add_assignment_rule('cmic', 'x_1 + x_2 + x_3')
#
model_net.add_species('doc', 'compartmentOne')
model_net.add_assignment_rule('doc', 'x_6 + x_8')
#
model_net.add_species('toc', 'compartmentOne')
model_net.add_assignment_rule('toc', 'x_1 + x_2 + x_3 + x_6 + x_7 + x_8 + x_9')

## Add inital activity fractions to parameters
model_net.add_parameter('rB0', 0.1)
model_net.add_parameter('rF0', 0.1)
#
model_net.set_var_ic('x_4', 'rB0')
model_net.set_var_ic('x_5', 'rF0')
#
model_net.add_assignment_rule('x_4', ' Phi21*x_8')
#
model_net.add_species('co2_tot', 'compartmentOne')
model_net.add_assignment_rule('co2_tot', 'x_10 + x_11*(1-YLHiq*(((cL + (t/(t**2 + bL))**3) - cL)/(cL + (1/2*sqrt(bL))**3)) - YLLoq*(1-(((cL + (t/(t**2 + bL))**3) - cL)/(cL + (1/2*sqrt(bL))**3))))*(cL + (t/(t**2 + bL))**3)')

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
IO.eqns_TeX_file(model_net, 'model36.tex')

## Create the model
m = Model([PECCAT_experiment.expt], [model_net])
p0= m.get_params().copy()

p36 = Utility.load('p35.geodesics.bpkl')
#

p36 = np.delete(p36,9)
for i, (key,value) in enumerate(p0.items()):
    p0[i] = p36[i]

## Set prior ranges from value/prior_range to value*prior_range

res = Residuals.PriorInLog('YLHiq_prior', 'YLHiq', np.log(p0[12]), np.log(np.sqrt(p0[12])))
m.AddResidual(res)
res = Residuals.PriorInLog('YLLoq_prior', 'YLLoq', np.log(p0[13]), np.log(np.sqrt(p0[13])))
m.AddResidual(res)
res = Residuals.PriorInLog('YrB_prior', 'YrB', np.log(p0[14]), np.log(np.sqrt(p0[14])))
m.AddResidual(res)
res = Residuals.PriorInLog('YrF_prior', 'YrF', np.log(p0[15]), np.log(np.sqrt(p0[15])))
m.AddResidual(res)
res = Residuals.PriorInLog('YRFP_prior', 'YRFP', np.log(p0[16]), np.log(np.sqrt(p0[16])))
m.AddResidual(res)
res = Residuals.PriorInLog('YsBHiq_prior', 'YsBHiq', np.log(p0[17]), np.log(np.sqrt(p0[17])))
m.AddResidual(res)#
res = Residuals.PriorInLog('YsBLoq_prior', 'YsBLoq', np.log(p0[18]), np.log(np.sqrt(p0[18])))
m.AddResidual(res)
res = Residuals.PriorInLog('YsBPHiq_prior', 'YsBPHiq', np.log(p0[19]), np.log(np.sqrt(p0[19])))
m.AddResidual(res)
res = Residuals.PriorInLog('YsBPLoq_prior', 'YsBPLoq', np.log(p0[20]), np.log(np.sqrt(p0[20])))
m.AddResidual(res)
res = Residuals.PriorInLog('YsBPP_prior', 'YsBPP', np.log(p0[21]), np.log(np.sqrt(p0[21])))
m.AddResidual(res)
res = Residuals.PriorInLog('YsFHiq_prior', 'YsFHiq', np.log(p0[22]), np.log(np.sqrt(p0[22])))
m.AddResidual(res)
res = Residuals.PriorInLog('YsFLoq_prior', 'YsFLoq', np.log(p0[23]), np.log(np.sqrt(p0[23])))
m.AddResidual(res)
res = Residuals.PriorInLog('rB0_prior', 'rB0', np.log(p0[24]), np.log(np.sqrt(p0[24])))
m.AddResidual(res)
res = Residuals.PriorInLog('rF0_prior', 'rF0', np.log(p0[25]), np.log(np.sqrt(p0[25])))
m.AddResidual(res)

## Optimize to fit data
print 'Initial Cost:', m.cost(p0)
popt = Optimization.fmin_lm_log_params(m, p0, maxiter=5, disp=True)
# Then we run Levenberg-Marquardt
#popt = Optimization.leastsq_log_params(m, popt1)
cost_opt = m.cost(popt)
print(popt)
print 'Optimized Cost:', m.cost(popt)

Utility.save(popt, 'popt36.model.bpkl')
