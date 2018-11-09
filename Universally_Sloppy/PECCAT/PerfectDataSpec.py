import scipy
from SloppyCell.ReactionNetworks import *
import os
from Nets import *


net.set_var_optimizable('Kf', False)
net.set_var_optimizable('nf', False)
net.set_var_optimizable('nc', False)
net.set_var_optimizable('theta', False)
net.set_var_optimizable('rhoB', False)
net.set_var_optimizable('cL', False)
net.set_var_optimizable('bL', False)


sens_traj = Dynamics.integrate_sensitivity(net, int_times, rtol=1e-9,fill_traj=True)

data_ids = [id for id in net.species.keys()
            if not net.get_var_constant(id)]

h, h_d = PerfectData.hessian_log_params(sens_traj, fixed_sf=True, return_dict=True, data_ids=data_ids)

scipy.savetxt('hessian.dat', h)
