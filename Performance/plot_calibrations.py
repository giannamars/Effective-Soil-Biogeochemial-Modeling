from SloppyCell.ReactionNetworks import *
import PECCAT_experiment
import PECCAT_experiment32_1
import PECCAT_experiment_MCPA
import numpy as np
import scipy

import matplotlib.pyplot as plt
import matplotlib

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
#
model_net.add_species('cLtot', 'cell')
model_net.add_assignment_rule('cLtot', '(cL + (t/(t**2 + bL))**3)*x_12')
#
model_net.add_species('cLhiq', 'cell')
model_net.add_assignment_rule('cLhiq', '(cL + (t/(t**2 + bL))**3)*x_12*YLHiq*(((cL + (t/(t**2 + bL))**3) - cL)/(cL + (1/2*sqrt(bL))**3))')
#
model_net.add_species('cLloq', 'cell')
model_net.add_assignment_rule('cLloq', '(cL + (t/(t**2 + bL))**3)*x_12*YLLoq*(1-(((cL + (t/(t**2 + bL))**3) - cL)/(cL + (1/2*sqrt(bL))**3)))')

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


## Output latex'ed equations for debugging
IO.eqns_TeX_file(model_net, 'model0.tex')

popt = Utility.load('popt00.model.bpkl')

m = Model([PECCAT_experiment.expt], [model_net])



## Read model32 from SBML file
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
#
model32_net.add_assignment_rule('x_4', ' Phi21*x_8')
#
model32_net.add_species('co2_tot', 'compartmentOne')
model32_net.add_assignment_rule('co2_tot', 'x_10 + x_11*(1-YLHiq*(((cL + (t/(t**2 + bL))**3) - cL)/(cL + (1/2*sqrt(bL))**3)) - YLLoq*(1-(((cL + (t/(t**2 + bL))**3) - cL)/(cL + (1/2*sqrt(bL))**3))))*(cL + (t/(t**2 + bL))**3)')

## Output latex'ed equations for debugging
IO.eqns_TeX_file(model32_net, 'model3232.tex')


popt32 = Utility.load('popt32.model.bpkl')
#
m32 = Model([PECCAT_experiment32_1.expt], [model32_net])


##
model32m_net = IO.from_SBML_file('PECCAT32m.xml', 'base')

## define bulk measurements
model32m_net.add_species('cmic', 'compartmentOne')
model32m_net.add_assignment_rule('cmic', 'x_1 + x_2 + x_3')
#
model32m_net.add_species('cBtot', 'compartmentOne')
model32m_net.add_assignment_rule('cBtot', 'x_1 + x_2')
#
model32m_net.add_species('doc', 'compartmentOne')
model32m_net.add_assignment_rule('doc', 'x_6 + x_8')
#
model32m_net.add_species('toc', 'compartmentOne')
model32m_net.add_assignment_rule('toc', 'x_1 + x_2 + x_3 + x_6 + x_7 + x_8 + x_9')


## Add inital activity fractions to parameters
model32m_net.add_parameter('rF0', 0.1)
#
model32m_net.set_var_ic('x_5', 'rF0')
#
model32m_net.add_assignment_rule('x_4', ' Phi21*x_8')




popt32mcpa = Utility.load('poptnew.model.bpkl')
m32mcpa = Model([PECCAT_experiment_MCPA.expt], [model32m_net])







# Plotting


fig, ax_traj = plt.subplots(figsize=fig_size, dpi=300)
fig.subplots_adjust(left=.18, bottom=.13, right=.99, top=.97)

traj = Dynamics.integrate(model_net, [0, 25],rtol=1e-9, params=popt, fill_traj=True)
ax_traj.plot(traj.get_times(), traj.get_var_traj('cBtot'), color='#525252', linewidth=2)
expt = m.exptColl.values()[0]
tt = expt.data['base']['cBtot'].keys()
yy = [expt.data['base']['cBtot'][t][0] for t in tt]
yerr = [expt.data['base']['cBtot'][t][1] for t in tt]
ax_traj.errorbar(tt,yy,yerr,fmt='.', elinewidth=1,color='xkcd:black', ms=7, label='MCPA+Litter')

traj32 = Dynamics.integrate(model32_net, [0, 25],rtol=1e-9, params=popt32, fill_traj=True)
ax_traj.plot(traj32.get_times(), traj32.get_var_traj('cBtot'), color='#e41a1c',  linestyle='--', linewidth=2)


traj32mcpa = Dynamics.integrate(model32m_net, [0, 25],rtol=1e-9, params=popt32mcpa, fill_traj=True)
ax_traj.plot(traj32mcpa.get_times(), traj32mcpa.get_var_traj('cBtot'), color='#4daf4a',  linestyle=':', linewidth=2)
expt = m32mcpa.exptColl.values()[0]
tt = expt.data['base']['cBtot'].keys()
yy = [expt.data['base']['cBtot'][t][0] for t in tt]
yerr = [expt.data['base']['cBtot'][t][1] for t in tt]
ax_traj.errorbar(tt,yy,yerr,fmt='s',ms=4, elinewidth=1,color='xkcd:black',markerfacecolor='none', markeredgecolor='xkcd:black', label='MCPA')

ax_traj.set_ylabel('Bacterial C [mg/g]')
#ax_traj.set_xlabel('Time [days]')
ax_traj.legend(frameon=False)
fig.savefig('cBtot32.pdf')
#plt.show()


fig, ax_traj = plt.subplots(figsize=fig_size, dpi=300)
#fig.subplots_adjust(left=.19, bottom=.16, right=.99, top=.9)
fig.subplots_adjust(left=.205, bottom=.13, right=.99, top=.97)


import matplotlib.ticker as ticker
#traj = Dynamics.integrate(model_net, [0, 25],rtol=1e-9, params=popt, fill_traj=True)
ax_traj.plot(traj.get_times(), traj.get_var_traj('x_2'), color='#525252', linewidth=2, label='Full (MCPA+Litter)')
expt = m.exptColl.values()[0]
tt = expt.data['base']['x_2'].keys()
yy = [expt.data['base']['x_2'][t][0] for t in tt]
yerr = [expt.data['base']['x_2'][t][1] for t in tt]
ax_traj.errorbar(tt,yy,yerr,fmt='.', elinewidth=1,color='xkcd:black', ms=7)
ax_traj.plot(traj32.get_times(), traj32.get_var_traj('x_2'), color='#e41a1c', linestyle='--', linewidth=2, label='Reduced (MCPA+Litter)')



traj32mcpa = Dynamics.integrate(model32m_net, [0, 25],rtol=1e-9, params=popt32mcpa, fill_traj=True)
ax_traj.plot(traj32mcpa.get_times(), traj32mcpa.get_var_traj('x_2'), color='#4daf4a',  linestyle=':', linewidth=2, label='Reduced (MCPA)')
expt = m32mcpa.exptColl.values()[0]
tt = expt.data['base']['x_2'].keys()
yy = [expt.data['base']['x_2'][t][0] for t in tt]
yerr = [expt.data['base']['x_2'][t][1] for t in tt]
ax_traj.errorbar(tt,yy,yerr,fmt='s', ms=4, elinewidth=1,color='xkcd:black',markerfacecolor='none', markeredgecolor='xkcd:black')

ax_traj.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0e'))
ax_traj.set_ylim(0.0e-4, 7e-4)
ax_traj.legend(frameon=False, loc='upper left')
ax_traj.set_ylabel('Spec. degrader C [mg/g]')
#ax_traj.set_xlabel('Time [days]')

fig.savefig('x232.pdf')





#fig, ax_traj = plt.subplots(figsize=(6.8,4.2), dpi=300)
fig, ax_traj = plt.subplots(figsize=fig_size, dpi=300)
#fig.subplots_adjust(left=.185, bottom=.16, right=.99, top=.97)
fig.subplots_adjust(left=.18, bottom=.13, right=.99, top=.97)


#traj = Dynamics.integrate(model_net, [0, 25],rtol=1e-9, params=popt, fill_traj=True)
ax_traj.plot(traj.get_times(), traj.get_var_traj('x_3'), color='#525252', linewidth=2)
expt = m.exptColl.values()[0]
tt = expt.data['base']['x_3'].keys()
yy = [expt.data['base']['x_3'][t][0] for t in tt]
yerr = [expt.data['base']['x_3'][t][1] for t in tt]
ax_traj.errorbar(tt,yy,yerr,fmt='.', elinewidth=1,color='xkcd:black', ms=7)
ax_traj.plot(traj32.get_times(), traj32.get_var_traj('x_3'), color='#e41a1c', linestyle='--', linewidth=2)



traj32mcpa = Dynamics.integrate(model32m_net, [0, 25],rtol=1e-9, params=popt32mcpa, fill_traj=True)
ax_traj.plot(traj32mcpa.get_times(), traj32mcpa.get_var_traj('x_3'), color='#4daf4a',  linestyle=':', linewidth=2)
expt = m32mcpa.exptColl.values()[0]
tt = expt.data['base']['x_3'].keys()
yy = [expt.data['base']['x_3'][t][0] for t in tt]
yerr = [expt.data['base']['x_3'][t][1] for t in tt]
ax_traj.errorbar(tt,yy,yerr,fmt='s', ms=4, elinewidth=1,color='xkcd:black',markerfacecolor='none', markeredgecolor='xkcd:black')


ax_traj.set_ylabel('Fungal C [mg/g]')
#ax_traj.set_xlabel('Time [days]')
fig.savefig('x332.pdf')

#fig, ax_traj = plt.subplots(figsize=(6.8,4.2), dpi=300)

fig, ax_traj = plt.subplots(figsize=fig_size, dpi=300)
#fig.subplots_adjust(left=.17, bottom=.16, right=.99, top=.97)
fig.subplots_adjust(left=.195, bottom=.13, right=.99, top=.97)


#traj = Dynamics.integrate(model_net, [0, 25],rtol=1e-9, params=popt, fill_traj=True)
ax_traj.plot(traj.get_times(), traj.get_var_traj('doc'), color='#525252', linewidth=2)
expt = m.exptColl.values()[0]
tt = expt.data['base']['doc'].keys()
yy = [expt.data['base']['doc'][t][0] for t in tt]
yerr = [expt.data['base']['doc'][t][1] for t in tt]
ax_traj.errorbar(tt,yy,yerr,fmt='.', elinewidth=1,color='xkcd:black', ms=7)
ax_traj.plot(traj32.get_times(), traj32.get_var_traj('doc'), color='#e41a1c', linestyle='--', linewidth=2)


traj32mcpa = Dynamics.integrate(model32m_net, [0, 25],rtol=1e-9, params=popt32mcpa, fill_traj=True)
ax_traj.plot(traj32mcpa.get_times(), traj32mcpa.get_var_traj('doc'), color='#4daf4a',  linestyle=':', linewidth=2)
expt = m32mcpa.exptColl.values()[0]
tt = expt.data['base']['doc'].keys()
yy = [expt.data['base']['doc'][t][0] for t in tt]
yerr = [expt.data['base']['doc'][t][1] for t in tt]
ax_traj.errorbar(tt,yy,yerr,fmt='s', ms=4, elinewidth=1,color='xkcd:black',markerfacecolor='none', markeredgecolor='xkcd:black')

ax_traj.set_ylabel('DOC [mg/g]')
#ax_traj.set_xlabel('Time [days]')

fig.savefig('doc32.pdf')

#fig, ax_traj = plt.subplots(figsize=(6.8,4.2), dpi=300)
fig, ax_traj = plt.subplots(figsize=fig_size, dpi=300)
#fig.subplots_adjust(left=.205, bottom=.16, right=.99, top=.97)
fig.subplots_adjust(left=.195, bottom=.13, right=.99, top=.97)


#traj = Dynamics.integrate(model_net, [0, 25],rtol=1e-9, params=popt, fill_traj=True)
ax_traj.plot(traj.get_times(), traj.get_var_traj('x_8'), color='#525252', linewidth=2)
expt = m.exptColl.values()[0]
tt = expt.data['base']['x_8'].keys()
yy = [expt.data['base']['x_8'][t][0] for t in tt]
yerr = [expt.data['base']['x_8'][t][1] for t in tt]
ax_traj.errorbar(tt,yy,yerr,fmt='.',ms=7, elinewidth=1,color='xkcd:black')
ax_traj.plot(traj32.get_times(), traj32.get_var_traj('x_7'), color='#e41a1c', linestyle='--', linewidth=2)


traj32mcpa = Dynamics.integrate(model32m_net, [0, 25],rtol=1e-9, params=popt32mcpa, fill_traj=True)
ax_traj.plot(traj32mcpa.get_times(), traj32mcpa.get_var_traj('x_7'), color='#4daf4a',  linestyle=':', linewidth=2)
expt = m32mcpa.exptColl.values()[0]
tt = expt.data['base']['x_7'].keys()
yy = [expt.data['base']['x_7'][t][0] for t in tt]
yerr = [expt.data['base']['x_7'][t][1] for t in tt]
ax_traj.errorbar(tt,yy,yerr,fmt='s', ms=4, elinewidth=1,color='xkcd:black',markerfacecolor='none', markeredgecolor='xkcd:black')

ax_traj.set_ylabel('MCPA-C [mg/g]')
#ax_traj.set_xlabel('Time [days]')

ax_traj.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0e'))

ax_traj.yaxis.labelpad = 0.1

fig.savefig('x832.pdf')

#fig, ax_traj = plt.subplots(figsize=(6.8,4.2), dpi=300)
fig, ax_traj = plt.subplots(figsize=fig_size, dpi=300)
#fig.subplots_adjust(left=.14, bottom=.16, right=.99, top=.97)
fig.subplots_adjust(left=.18, bottom=.13, right=.99, top=.97)


#traj = Dynamics.integrate(model_net, [0, 25],rtol=1e-9, params=popt, fill_traj=True)
ax_traj.plot(traj.get_times(), traj.get_var_traj('toc'), color='#525252', linewidth=2)
expt = m.exptColl.values()[0]
tt = expt.data['base']['toc'].keys()
yy = [expt.data['base']['toc'][t][0] for t in tt]
yerr = [expt.data['base']['toc'][t][1] for t in tt]
ax_traj.errorbar(tt,yy,yerr,fmt='.',ms=7, elinewidth=1,color='xkcd:black')
ax_traj.plot(traj32.get_times(), traj32.get_var_traj('toc'), color='#e41a1c', linestyle='--', linewidth=2)

traj32mcpa = Dynamics.integrate(model32m_net, [0, 25],rtol=1e-9, params=popt32mcpa, fill_traj=True)
ax_traj.plot(traj32mcpa.get_times(), traj32mcpa.get_var_traj('toc'), color='#4daf4a',  linestyle=':', linewidth=2)
expt = m32mcpa.exptColl.values()[0]
tt = expt.data['base']['toc'].keys()
yy = [expt.data['base']['toc'][t][0] for t in tt]
yerr = [expt.data['base']['toc'][t][1] for t in tt]
ax_traj.errorbar(tt,yy,yerr,fmt='s', ms=4, elinewidth=1,color='xkcd:black',markerfacecolor='none', markeredgecolor='xkcd:black')


ax_traj.set_ylabel('TOC [mg/g]')
#ax_traj.set_xlabel('Time [days]')

fig.savefig('toc32.pdf')

#fig, ax_traj = plt.subplots(figsize=(6.8,4.2), dpi=300)
fig, ax_traj = plt.subplots(figsize=fig_size, dpi=300)
#fig.subplots_adjust(left=.17, bottom=.16, right=.99, top=.97)
fig.subplots_adjust(left=.18, bottom=.13, right=.99, top=.97)

#traj = Dynamics.integrate(model_net, [0, 25],rtol=1e-9, params=popt, fill_traj=True)
ax_traj.plot(traj.get_times(), traj.get_var_traj('cmic'), color='#525252', linewidth=2)
expt = m.exptColl.values()[0]
tt = expt.data['base']['cmic'].keys()
yy = [expt.data['base']['cmic'][t][0] for t in tt]
yerr = [expt.data['base']['cmic'][t][1] for t in tt]
ax_traj.errorbar(tt,yy,yerr,fmt='.',ms=7, elinewidth=1,color='xkcd:black')
ax_traj.plot(traj32.get_times(), traj32.get_var_traj('cmic'), color='#e41a1c', linestyle='--', linewidth=2)


traj32mcpa = Dynamics.integrate(model32m_net, [0, 25],rtol=1e-9, params=popt32mcpa, fill_traj=True)
ax_traj.plot(traj32mcpa.get_times(), traj32mcpa.get_var_traj('cmic'), color='#4daf4a',  linestyle=':', linewidth=2)
expt = m32mcpa.exptColl.values()[0]
tt = expt.data['base']['cmic'].keys()
yy = [expt.data['base']['cmic'][t][0] for t in tt]
yerr = [expt.data['base']['cmic'][t][1] for t in tt]
ax_traj.errorbar(tt,yy,yerr,fmt='s', ms=4, elinewidth=1,color='xkcd:black',markerfacecolor='none', markeredgecolor='xkcd:black')


ax_traj.set_ylabel('Microbial C [mg/g]')
#ax_traj.set_xlabel('Time [days]')
fig.savefig('cmic32.pdf')

#fig, ax_traj = plt.subplots(figsize=(6.8,4.2), dpi=300)
fig, ax_traj = plt.subplots(figsize=fig_size, dpi=300)
#fig.subplots_adjust(left=.165, bottom=.16, right=.99, top=.97)
fig.subplots_adjust(left=.165, bottom=.13, right=.99, top=.97)

#traj = Dynamics.integrate(model_net, [0, 25],rtol=1e-9, params=popt, fill_traj=True)
ax_traj.plot(traj.get_times(), traj.get_var_traj('co2_tot'), color='#525252', linewidth=2)
expt = m.exptColl.values()[0]
tt = expt.data['base']['co2_tot'].keys()
yy = [expt.data['base']['co2_tot'][t][0] for t in tt]
yerr = [expt.data['base']['co2_tot'][t][1] for t in tt]
ax_traj.errorbar(tt,yy,yerr,fmt='.',ms=7, elinewidth=1,color='xkcd:black')
ax_traj.plot(traj32.get_times(), traj32.get_var_traj('co2_tot'), color='#e41a1c', linestyle='--', linewidth=2)



traj32mcpa = Dynamics.integrate(model32m_net, [0, 25],rtol=1e-9, params=popt32mcpa, fill_traj=True)
ax_traj.plot(traj32mcpa.get_times(), traj32mcpa.get_var_traj('x_10'), color='#4daf4a',  linestyle=':', linewidth=2)
expt = m32mcpa.exptColl.values()[0]
tt = expt.data['base']['x_10'].keys()
yy = [expt.data['base']['x_10'][t][0] for t in tt]
yerr = [expt.data['base']['x_10'][t][1] for t in tt]
ax_traj.errorbar(tt,yy,yerr,fmt='s', ms=4, elinewidth=1,color='xkcd:black',markerfacecolor='none', markeredgecolor='xkcd:black')

ax_traj.set_ylabel(r'$CO_2$-C [mg/g]')
#ax_traj.set_xlabel('Time [days]')

fig.savefig('co2tot32.pdf')


#fig, ax_traj = plt.subplots(figsize=(6.8,4.2), dpi=300)

#traj = Dynamics.integrate(model_net, [0, 25],rtol=1e-9, params=popt, fill_traj=True)
#ax_traj.plot(traj.get_times(), traj.get_var_traj('cLtot'), color='xkcd:black', linewidth=2)
#ax_traj.plot(traj.get_times(), traj.get_var_traj('cLhiq'), color='xkcd:black', linewidth=2, linestyle='--')
#ax_traj.plot(traj.get_times(), traj.get_var_traj('cLloq'), color='xkcd:black', linewidth=2, linestyle=':')
#fig.savefig('cL.svg')



plt.show()
