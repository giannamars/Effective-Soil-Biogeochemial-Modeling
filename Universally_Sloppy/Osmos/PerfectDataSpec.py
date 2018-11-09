import scipy
from SloppyCell.ReactionNetworks import *
import os
from Nets import *

net.set_var_optimizable('b', False)
net.set_var_optimizable('w', False)



for id in net.optimizableVars.keys():
    if net.get_var_ic(id)==0:
        net.set_var_optimizable(id,False)



sens_traj = Dynamics.integrate_sensitivity(net, int_times, rtol=1e-9,fill_traj=True)

data_ids = [id for id in net.species.keys()
            if not net.get_var_constant(id)]




h, h_d = PerfectData.hessian_log_params(sens_traj, fixed_sf=True, return_dict=True, data_ids=data_ids)

scipy.savetxt('hessian.dat', h)
