from SloppyCell.ReactionNetworks import *
import numpy as np
import scipy
import german_expt

model_net = IO.from_SBML_file('germans1.xml', 'base')
model_net.set_var_optimizable('I', False)

# Define respiration
model_net.add_species('Resp', 'compartmentOne')
model_net.add_assignment_rule('Resp', '1 - eps*Phi1*x_1*x_2')

## Output latex'ed equations for debugging
IO.eqns_TeX_file(model_net, 'model1.tex')


# Create the model
m = Model([german_expt.expt], [model_net])
p0= m.get_params().copy()
print(len(p0))


## Set prior ranges from value/prior_range to value*prior_range in interval [0,1] for substrate utilization efficiencies
res = Residuals.PriorInLog('eps_prior', 'eps', np.log(p0[0]), np.log(np.sqrt(p0[0])))
m.AddResidual(res)

#p1 = Utility.load('p0.geodesics.bpkl')
#p1[2] = p1[3]/p1[1]
#p1 = np.delete(p1, 3)
#
#
#for i, (key,value) in enumerate(p0.items()):
#    p0[i] = p1[i]
p0[1] = 3.5e-5

## Optimize to fit data
print 'Initial Cost:', m.cost(p0)
popt = Optimization.fmin_lm_log_params(m, p0, maxiter=100, disp=True)
cost_opt = m.cost(popt)
print(popt)
print 'Optimized Cost:', m.cost(popt)
Utility.save(popt, 'popt1.model.bpkl')



#import matplotlib.pyplot as plt
#import matplotlib
#
#fig_width = 3  # inches
#golden_mean = (np.sqrt(5)-1.0)/2.0   # aethestic ratio
#fig_height = fig_width*golden_mean
#dpi = 300.0
#
#fig_size = [fig_width, fig_height]
#
#params = {'backend': 'ps',
#          'axes.labelsize': 8,
#          'axes.titlesize': 8,
#          'font.size': 8,
#          'legend.fontsize': 8,
#          'xtick.labelsize': 8,
#          'ytick.labelsize': 8,
#          'text.usetex': False,
#          'figure.figsize': fig_size}
#
#
#plt.rcParams.update(params)
#
#fig, ax_traj = plt.subplots(figsize=fig_size, dpi=300)
#fig.subplots_adjust(left=.18, bottom=.13, right=.99, top=.97)
#
#traj = Dynamics.integrate(model_net, [0, 8640],rtol=1e-9, params=popt, fill_traj=True)
#ax_traj.plot(traj.get_times(), traj.get_var_traj('Resp'), color='red', linewidth=2)
#expt = m.exptColl.values()[0]
#tt = expt.data['base']['Resp'].keys()
#yy = [expt.data['base']['Resp'][t][0] for t in tt]
#yerr = [expt.data['base']['Resp'][t][1] for t in tt]
#ax_traj.errorbar(tt,yy,yerr,fmt='.', elinewidth=1,color='xkcd:black', ms=7, label='Perfect')
#
#ax_traj.set_ylabel('Resp')
#ax_traj.set_xlabel('Time')
#fig.savefig('german_Resp_0_fit.pdf')


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
