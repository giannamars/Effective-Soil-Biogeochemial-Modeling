from SloppyCell.ReactionNetworks import *
import numpy as np
import scipy
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.ticker as ticker


fig_width = 3  # inches
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


plt.rcParams.update(params)



## Read  model 0 from SBML file
model_net = IO.from_SBML_file('PECCAT0.xml', 'base')
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
model_net.add_parameter('rB0', 0.1)
model_net.add_parameter('rBP0', 0.2)
model_net.add_parameter('rF0', 0.1)
#
model_net.set_var_ic('x_4', 'rB0')
model_net.set_var_ic('x_5', 'rBP0')
model_net.set_var_ic('x_5', 'rF0')
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

# Perturb Model 0
modelp_net = model_net.copy('perturbed')
modelp_net.set_var_ic('x_12', 0.0)

# Calculate Model36 perturbed predictions

tt = np.linspace(0,25,100)

def prune_ens(ens, Npt):
   return ens[1::len(ens)//Npt]
#    return ens[1::10000]

def calc_trajectories(ens,tt, net):
    yys = []
    for p in ens:
	traj = Dynamics.integrate(net, tt, rtol=1e-9, params=p, fill_traj=False)
        yy = traj.get_var_traj('x_8')
        yys.append(yy)
    return yys

def calc_trajectory_bounds(yys, lw_percent=2.5, up_percent=97.5):
    up, lw = [], []
    for yt in np.array(yys).transpose():
        up.append(scipy.stats.scoreatpercentile(yt, 97.5))
        lw.append(scipy.stats.scoreatpercentile(yt, 2.5))
    return np.array(lw), np.array(up)

# calculate ensemble predictions for perturbed net
ens_data = Utility.load('model0.ens_data_mio.bpkl')


print(ens_data.shape)
ens_pruned = prune_ens(ens_data, 1000)
print(ens_pruned.shape)
print(ens_pruned[0].shape)
pred_data = calc_trajectories(ens_pruned, tt, modelp_net)
pred_data_lw, pred_data_up = calc_trajectory_bounds(pred_data)



## Read  model 32 from SBML file
model32_net = IO.from_SBML_file('PECCAT32.xml', 'base')
## define bulk measurements
model32_net.add_species('cmic', 'compartmentOne')
model32_net.add_assignment_rule('cmic', 'x_1 + x_2 + x_3')
#
model32_net.add_species('cBtot', 'compartmentOne')
model32_net.add_assignment_rule('cBtot', 'x_1 + x_2')
#
model32_net.add_species('doc', 'compartmentOne')
model32_net.add_assignment_rule('doc', 'x_6 + x_8')
#
model32_net.add_species('toc', 'compartmentOne')
model32_net.add_assignment_rule('toc', 'x_1 + x_2 + x_3 + x_6 + x_7 + x_8 + x_9')

## Add inital activity fractions to parameters
model32_net.add_parameter('rF0', 0.1)
#
model32_net.set_var_ic('x_5', 'rF0')
model32_net.add_assignment_rule('x_4', 'Phi21*x_8')
#
model32_net.add_species('co2_tot', 'compartmentOne')
model32_net.add_assignment_rule('co2_tot', 'x_10 + x_11*(1-YLHiq*(((cL + (t/(t**2 + bL))**3) - cL)/(cL + (1/2*sqrt(bL))**3)) - YLLoq*(1-(((cL + (t/(t**2 + bL))**3) - cL)/(cL + (1/2*sqrt(bL))**3))))*(cL + (t/(t**2 + bL))**3)')

## Don't optimize MCPA sorption kinetics and soil properties
model32_net.set_var_optimizable('Kf', False)
model32_net.set_var_optimizable('nf', False)
model32_net.set_var_optimizable('nc', False)
model32_net.set_var_optimizable('theta', False)
model32_net.set_var_optimizable('rhoB', False)
model32_net.set_var_optimizable('cL', False)
model32_net.set_var_optimizable('bL', False)
model32_net.set_var_optimizable('YLHiq', True)
model32_net.set_var_optimizable('YLLoq', True)

# Perturb Model 0
modelp32_net = model32_net.copy('perturbed32')
modelp32_net.set_var_ic('x_11', 0.0)

# Calculate Model36 perturbed predictions

def calc_trajectories32(ens,tt, net):
    yys = []
    for p in ens:
        traj = Dynamics.integrate(net, tt, rtol=1e-9, params=p, fill_traj=False)
        yy = traj.get_var_traj('x_7')
        yys.append(yy)
    return yys


# calculate ensemble predictions for perturbed net
ens_data32 = Utility.load('model32.ens_data1.bpkl')

print(ens_data32.shape)
ens_pruned32 = prune_ens(ens_data32, 10)
print(ens_pruned32.shape)
print(ens_pruned32[0].shape)
pred_data32 = calc_trajectories32(ens_pruned32, tt, modelp32_net)
pred_data_lw32, pred_data_up32 = calc_trajectory_bounds(pred_data32)


# Calulate Model 0 reverse predictions

model_net.set_var_ic('x_12', '0')
# Perturb Model 0
modelp0_net = model_net.copy('perturbed')
modelp0_net.set_var_ic('x_12', 0.089402100177517582)

# Calculate Model00  perturbed predictions

tt = np.linspace(0,25,100)

def prune_ens(ens, Npt):
   return ens[1::len(ens)//Npt]
#    return ens[1::10000]

def calc_trajectories(ens,tt, net):
    yys = []
    for p in ens:
	traj = Dynamics.integrate(net, tt, rtol=1e-9, params=p, fill_traj=False)
        yy = traj.get_var_traj('x_8')
        yys.append(yy)
    return yys

def calc_trajectory_bounds(yys, lw_percent=2.5, up_percent=97.5):
    up, lw = [], []
    for yt in np.array(yys).transpose():
        up.append(scipy.stats.scoreatpercentile(yt, 97.5))
        lw.append(scipy.stats.scoreatpercentile(yt, 2.5))
    return np.array(lw), np.array(up)

# calculate ensemble predictions for perturbed net
ens_data0 = Utility.load('model00.ens_data.bpkl')


print(ens_data0.shape)
ens_pruned0 = prune_ens(ens_data0, 10)
print(ens_pruned0.shape)
print(ens_pruned0[0].shape)
pred_data0 = calc_trajectories(ens_pruned0, tt, modelp0_net)
pred_data_lw0, pred_data_up0 = calc_trajectory_bounds(pred_data0)


## Calculate Model32 reverse prediction

model32_net.set_var_ic('x_11', '0')
# Perturb Model 320
modelp320_net = model32_net.copy('perturbed')
modelp320_net.set_var_ic('x_11', 0.089402100177517582)

# Calculate Model320  perturbed predictions

tt = np.linspace(0,25,100)

def prune_ens(ens, Npt):
   return ens[1::len(ens)//Npt]
#   return ens[1::1000]
def calc_trajectories(ens,tt, net):
    yys = []
    for p in ens:
	traj = Dynamics.integrate(net, tt, rtol=1e-9, params=p, fill_traj=False)
        yy = traj.get_var_traj('x_7')
        yys.append(yy)
    return yys

def calc_trajectory_bounds(yys, lw_percent=2.5, up_percent=97.5):
    up, lw = [], []
    for yt in np.array(yys).transpose():
        up.append(scipy.stats.scoreatpercentile(yt, 97.5))
        lw.append(scipy.stats.scoreatpercentile(yt, 2.5))
    return np.array(lw), np.array(up)

# calculate ensemble predictions for perturbed net
ens_data320 = Utility.load('model320.ens_data1.bpkl')


print(ens_data320.shape)
ens_pruned320 = prune_ens(ens_data320, 10)
print(ens_pruned320.shape)
print(ens_pruned320[0].shape)
pred_data320 = calc_trajectories(ens_pruned320, tt, modelp320_net)
pred_data_lw320, pred_data_up320 = calc_trajectory_bounds(pred_data320)




## Load data sets

import PECCAT_experiment
import PECCAT_experiment_MCPA

m = Model([PECCAT_experiment.expt], [model_net])
m32 = Model([PECCAT_experiment_MCPA.expt], [model32_net])





#Plotting

fig, ax = plt.subplots(figsize=fig_size, dpi=300)
fig.subplots_adjust(left=.195, bottom=.13, right=.99, top=.97)


ax.plot(tt, pred_data_lw, color='#525252', linewidth=2, label='Full (MCPA)')
ax.plot(tt, pred_data_up, color='#525252', linewidth=2)
#ax.fill_between(tt,pred_data_lw, pred_data_up, color='black')



ax.plot(tt, pred_data_lw32, color='#525252', ls='--', linewidth=2, label='Reduced (MCPA)')
ax.plot(tt, pred_data_up32, color='#525252', ls ='--', linewidth=2)
#ax.fill_between(tt,pred_data_lw32, pred_data_up32, color='black')
#ax.set_ylim(0,0.00002)




ax.plot(tt, pred_data_lw0, color='#8da0cb', ls='-', linewidth=2, label='Full (MCPA+Litter)')
ax.plot(tt, pred_data_up0, color='#8da0cb', ls='-',linewidth=2)
#ax.fill_between(tt,pred_data_lw, pred_data_up, color='black')



ax.plot(tt, pred_data_lw320, color='#8da0cb',ls='--', linewidth=2, label='Reduced (MCPA+Litter)')
ax.plot(tt, pred_data_up320, color='#8da0cb', ls='--', linewidth=2)
#ax.fill_between(tt,pred_data_lw, pred_data_up, color='black')


expt = m.exptColl.values()[0]
tt = expt.data['base']['x_8'].keys()
yy = [expt.data['base']['x_8'][t][0] for t in tt]
yerr = [expt.data['base']['x_8'][t][1] for t in tt]
ax.errorbar(tt,yy,yerr,fmt='.',ms=7, elinewidth=1,color='xkcd:black',)


expt = m32.exptColl.values()[0]
tt = expt.data['base']['x_7'].keys()
yy = [expt.data['base']['x_7'][t][0] for t in tt]
yerr = [expt.data['base']['x_7'][t][1] for t in tt]
ax.errorbar(tt,yy,yerr,fmt='s',ms=4, elinewidth=1,color='xkcd:black', markerfacecolor='none', markeredgecolor='black')


handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc='center right', frameon=False, fontsize=7)


ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0e'))

#ax.set_xlabel('Time [days]')
ax.set_ylabel('MCPA-C [mg/g]')
ax.yaxis.labelpad = 0.1
fig.savefig('model_predictions.pdf')





plt.show()
