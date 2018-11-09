from SloppyCell.ReactionNetworks import *
import PECCAT_experiment
import numpy as np
import scipy
import matplotlib.pyplot as plt

## Read model from SBML file
model_net = IO.from_SBML_file('PECCAT0.xml', 'base')

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
#model_net.add_parameter('Lit0', 0.1)
#model_net.set_var_ic('x_12', 'Lit0')
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
IO.eqns_TeX_file(model_net, 'model0.tex')

## Create the model
m = Model([PECCAT_experiment.expt], [model_net])
p0= m.get_params().copy()
print(len(p0))
p0 = Utility.load('popt0.model.bpkl')
#p0[0]=1.2
#p0[1]=1.2
#p0[2]=1.2
#p0[10]=30.0
#p0[11]=30.0
## Set prior ranges from value/prior_range to value*prior_range

prior_range = 100

for key, value in p0.items()[0:len(p0)-15]:

    res = Residuals.PriorInLog('%s_prior' % key, key, np.log(value), np.log(np.sqrt(prior_range)))
    m.AddResidual(res)


res = Residuals.PriorInLog('YLHiq_prior', 'YLHiq', np.log(p0[44]), np.log(np.sqrt(p0[44])))
m.AddResidual(res)
res = Residuals.PriorInLog('YLLoq_prior', 'YLLoq', np.log(p0[45]), np.log(np.sqrt(p0[45])))
m.AddResidual(res)
res = Residuals.PriorInLog('YrB_prior', 'YrB', np.log(p0[46]), np.log(np.sqrt(p0[46])))
m.AddResidual(res)
res = Residuals.PriorInLog('YrF_prior', 'YrF', np.log(p0[47]), np.log(np.sqrt(p0[47])))
m.AddResidual(res)
res = Residuals.PriorInLog('YRFP_prior', 'YRFP', np.log(p0[48]), np.log(np.sqrt(p0[48])))
m.AddResidual(res)
res = Residuals.PriorInLog('YsBHiq_prior', 'YsBHiq', np.log(p0[49]), np.log(np.sqrt(p0[49])))
m.AddResidual(res)#
res = Residuals.PriorInLog('YsBLoq_prior', 'YsBLoq', np.log(p0[50]), np.log(np.sqrt(p0[50])))
m.AddResidual(res)
res = Residuals.PriorInLog('YsBPHiq_prior', 'YsBPHiq', np.log(p0[51]), np.log(np.sqrt(p0[51])))
m.AddResidual(res)
res = Residuals.PriorInLog('YsBPLoq_prior', 'YsBPLoq', np.log(p0[52]), np.log(np.sqrt(p0[52])))
m.AddResidual(res)
res = Residuals.PriorInLog('YsBPP_prior', 'YsBPP', np.log(p0[53]), np.log(np.sqrt(p0[53])))
m.AddResidual(res)
res = Residuals.PriorInLog('YsFHiq_prior', 'YsFHiq', np.log(p0[54]), np.log(np.sqrt(p0[54])))
m.AddResidual(res)
res = Residuals.PriorInLog('YsFLoq_prior', 'YsFLoq', np.log(p0[55]), np.log(np.sqrt(p0[55])))
m.AddResidual(res)
res = Residuals.PriorInLog('rB0_prior', 'rB0', np.log(p0[56]), np.log(np.sqrt(p0[56])))
m.AddResidual(res)
res = Residuals.PriorInLog('rBP0_prior', 'rBP0', np.log(p0[57]), np.log(np.sqrt(p0[57])))
m.AddResidual(res)
res = Residuals.PriorInLog('rF0_prior', 'rF0', np.log(p0[58]), np.log(np.sqrt(p0[58])))
m.AddResidual(res)

## Optimize to fit data
print 'Initial Cost:', m.cost(p0)
popt = Optimization.fmin_lm_log_params(m, p0, maxiter=10, disp=True)
# Then we run Levenberg-Marquardt
#popt = Optimization.leastsq_log_params(m, popt1)
cost_opt = m.cost(popt)
print(popt)
print(len(popt))
print 'Optimized Cost:', m.cost(popt)

Utility.save(popt, 'popt00.model.bpkl')



## Ensemble fit
print('Calculating Hessian')
j = m.jacobian_log_params_sens(np.log(popt))
jtj = np.dot(np.transpose(j),j)
np.savetxt('hessian0.dat', jtj)
e,v = Utility.eig(jtj)
e = scipy.real(e)
#Plotting.figure(0)
#l = Plotting.plot_eigval_spectrum(e, offset=0.15, widths=0.7, lw=2)
#
print('Generating ensemble')
ens_data, costs, ratio = Ensembles.ensemble_log_params(m, popt, jtj, steps=30000)
#


ens_PCA = Ensembles.PCA_eig_log_params(ens_data)
eigsp = 1/ens_PCA[0]**2
np.savetxt('hessian.dat', eigsp)
#e,v = Utility.eig(ens_PCA[0])
#e = scipy.real(e)

Plotting.figure(1)
l = Plotting.plot_eigval_spectrum(eigsp)
plt.savefig('pcaspec.png')


#ens_data = np.array(ens_data)
#Utility.save(ens_data, 'model0.ens_data.bpkl')
#costs = np.array(costs)
#Utility.save(costs, 'model0.costs.bpkl')



## Plotting

#Plotting.figure(1)
#Plotting.plot_model_results(m, data_to_plot = 'cBtot')
#Plotting.figure(2)
##Plotting.plot_model_results(m, data_to_plot = 'x_2')
#Plotting.figure(3)
#Plotting.plot_model_results(m, data_to_plot = 'x_3')
#Plotting.figure(4)
#Plotting.plot_model_results(m, data_to_plot = 'doc')
#Plotting.figure(5)
#Plotting.plot_model_results(m, data_to_plot = 'x_8')
#Plotting.figure(6)
#Plotting.plot_model_results(m, data_to_plot = 'toc')
#Plotting.figure(7)
#Plotting.plot_model_results(m, data_to_plot = 'co2_tot')
#Plotting.figure(8)
#Plotting.plot_model_results(m, data_to_plot = 'cmic')



#traj = Dynamics.integrate(model_net, [0, 25],rtol=1e-9, params=popt, fill_traj=True)

#Plotting.figure(9)
#plt.plot(traj.get_times(), traj.get_var_traj('x_1'))

#Plotting.figure(10)
#plt.plot(traj.get_times(), traj.get_var_traj('x_2'))

#Plotting.figure(11)
#plt.plot(traj.get_times(), traj.get_var_traj('x_3'))

#Plotting.figure(12)
#plt.plot(traj.get_times(), traj.get_var_traj('x_4'))

#Plotting.figure(13)
#plt.plot(traj.get_times(), traj.get_var_traj('x_5'))


#Plotting.figure(14)
#plt.plot(traj.get_times(), traj.get_var_traj('x_6'))


#Plotting.figure(15)
#plt.plot(traj.get_times(), traj.get_var_traj('x_7'))


#Plotting.figure(16)
#plt.plot(traj.get_times(), traj.get_var_traj('x_8'))

#Plotting.figure(17)
#plt.plot(traj.get_times(), traj.get_var_traj('x_9'))

#Plotting.figure(18)
#plt.plot(traj.get_times(), traj.get_var_traj('x_10'))

#Plotting.figure(19)
#plt.plot(traj.get_times(), traj.get_var_traj('x_11'))

#Plotting.figure(20)
#plt.plot(traj.get_times(), traj.get_var_traj('x_12'))



plt.show()
