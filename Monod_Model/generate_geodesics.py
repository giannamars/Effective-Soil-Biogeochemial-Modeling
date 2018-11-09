from SloppyCell.ReactionNetworks import *
from numpy import *
import geodesics
import monod_model

x = log(monod_model.popt)
v = array([0,1])
tmax = 10
Avv = None

func = lambda logp: array(monod_model.m.res_log_params(logp))
jacobian = lambda logp: array(monod_model.m.jacobian_log_params_sens(logp).values())

u, vects = Utility.eig(monod_model.jtj)
# Draw orthogonal geodesics coming out from optimal parameter values in most sloppy direction.
theta0 = arctan2(vects[1,-1], vects[0,-1])
theta_list =  theta0 + linspace(0,2*pi,2,endpoint=False)

result_list = []
for theta in theta_list:
    v = array([cos(theta), sin(theta)])
    xs, vs, ts = geodesics.geodesic(x, v, tmax, func, jacobian, Avv,
                                    maxsteps=int(2e4), rtol=1e-6, atol=1e-6)
    result_list.append((xs,vs,ts))
    print('Got to t={0}'.format(ts[-1]))

Utility.save(result_list, 'monod.geodesics.bpkl')
