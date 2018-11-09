import scipy
from SloppyCell.ReactionNetworks import *
import os
from Nets import *

net.set_var_optimizable('I', False)


sens_traj = Dynamics.integrate_sensitivity(net, int_times, rtol=1e-9,fill_traj=True)

data_ids = [id for id in net.species.keys()
            if not net.get_var_constant(id)]

h, h_d = PerfectData.hessian_log_params(sens_traj, fixed_sf=True, return_dict=True, data_ids=data_ids)

scipy.savetxt('hessian.dat', h)

