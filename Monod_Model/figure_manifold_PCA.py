from numpy import *
from SloppyCell.ReactionNetworks import *
import monod_net, monod_expt
import monod_model


sample_times = linspace(0,40,6)

# Calculate a model output surface

Npt = 151
land_xx = logspace(-3, 2, Npt)
land_yy = logspace(-3, 2, Npt)
sample_data = empty((len(land_yy), len(land_xx), 2*(len(sample_times)-1)))

for ii,y in enumerate(land_yy):
    for jj,x in enumerate(land_xx):
        p = [x, y]

        traj = Dynamics.integrate(monod_net.net, sample_times, params=p, fill_traj=False)

        for kk, var in enumerate(['x_1', 'x_2']):
            start = (len(sample_times)-1)*kk
            end = (len(sample_times)-1)*(kk+1)
            sample_data[ii,jj,start:end] = abs(traj.get_var_traj(var)[1:])


print(sample_data.shape)

# Take PCA of sample_data, reshaped (101*101, 2*sample_times)

reshaped = sample_data.reshape((-1, sample_data.shape[-1]))
#print(reshaped.size)
#reshaped = reshaped[logical_not(any(isnan(reshaped), axis=-1))]
# u has 10 eigenvalues, v has columns of eigenvectors (10,10)

u, v = Ensembles.PCA_eig(reshaped)

savetxt('monod.model_manifold.nolog.PCA.us.dat', u)

# Project trajectories onto PCA axes
mean_sample_data = mean(reshaped, axis=0)
all_projs = dot(sample_data - mean_sample_data, v)
print(all_projs.shape)

## Projection of optimal trajectory
traj = Dynamics.integrate(monod_net.net, sample_times, params = monod_model.popt, fill_traj=False)
popt_sample_data = zeros(all_projs.shape[-1])
#print(popt_sample_data.shape)

for kk, var in enumerate(['x_1', 'x_2']):
    start = (len(sample_times)-1)*kk
    end = (len(sample_times)-1)*(kk+1)
    popt_sample_data[start:end] = abs(traj.get_var_traj(var)[1:])

popt_proj = dot(popt_sample_data - mean_sample_data, v)


## Projection of Data point
data_sample_data = zeros(all_projs.shape[-1])
for kk, var in enumerate(['x_1', 'x_2']):
    data = [val for (time, (val,sigma)) in sorted(monod_expt.expt.get_data()['monod'][var].items())]
    start = (len(sample_times)-1)*kk
    end = (len(sample_times)-1)*(kk+1)
    data_sample_data[start:end] = data

data_proj = dot(data_sample_data - mean_sample_data, v)



# Projection of geodesics
geodesics = Utility.load('monod.geodesics.bpkl')
all_gprojs = []
for xs,vs,ts in geodesics:
    gsamps = []
    for logp in xs:
        p = exp(logp)
        if not(1e-3 <= p[0] <= 1e2) or not (1e-3 <= p[1] <= 1e2):
            continue

        traj = Dynamics.integrate(monod_net.net, sample_times, params = p, fill_traj=False)
        sample_data = []
        for kk, var in enumerate(['x_1', 'x_2']):
            start = (len(sample_times)-1)*kk
            end = (len(sample_times)-1)*(kk+1)
            sample_data.extend(abs(traj.get_var_traj(var)[1:]))
        gsamps.append(sample_data)
    gsamps = array(gsamps)
    gprojs = dot(gsamps - mean_sample_data, v)
    all_gprojs.append(gprojs)



## Projection of p=[20,20] trajectory
traj = Dynamics.integrate(monod_net.net, sample_times, params = [2.5,0.2], fill_traj=False)
p20_sample_data = zeros(all_projs.shape[-1])
#print(popt_sample_data.shape)

for kk, var in enumerate(['x_1', 'x_2']):
    start = (len(sample_times)-1)*kk
    end = (len(sample_times)-1)*(kk+1)
    p20_sample_data[start:end] = abs(traj.get_var_traj(var)[1:])

p20_proj = dot(p20_sample_data - mean_sample_data, v)



## Projection of p=[1e-2,0.3*1e-1] trajectory
traj = Dynamics.integrate(monod_net.net, sample_times, params = [5,1.1], fill_traj=False)
p21_sample_data = zeros(all_projs.shape[-1])

for kk, var in enumerate(['x_1', 'x_2']):
    start = (len(sample_times)-1)*kk
    end = (len(sample_times)-1)*(kk+1)
    p21_sample_data[start:end] = abs(traj.get_var_traj(var)[1:])

p21_proj = dot(p21_sample_data - mean_sample_data, v)



## Projection of p=[1e-2,0.3*1e-1] trajectory
traj = Dynamics.integrate(monod_net.net, sample_times, params = [0.15,0.3], fill_traj=False)
p22_sample_data = zeros(all_projs.shape[-1])

for kk, var in enumerate(['x_1', 'x_2']):
    start = (len(sample_times)-1)*kk
    end = (len(sample_times)-1)*(kk+1)
    p22_sample_data[start:end] = abs(traj.get_var_traj(var)[1:])

p22_proj = dot(p22_sample_data - mean_sample_data, v)




## Plotting
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


plt.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
plt.rc('text', usetex=False)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', labelsize=8)
plt.rc('legend', fontsize=8)
width = 3.0
golden_mean = (sqrt(5)-1.0)/2.0   # aethestic ratio
height = width*golden_mean



fig = plt.figure()
fig.subplots_adjust(left=.12, bottom=.15, right=.99, top=.95)

ax3d = fig.gca(projection='3d')

land_xx, land_yy, Z = Utility.load('monod.model_surface.bpkl')

# Calculate cost as distance to data point
dist = empty((all_projs.shape[0], all_projs.shape[1]))
b = array((data_proj[-3], data_proj[-2], data_proj[-1] ))
for ii in range(all_projs.shape[0]):
    for jj in range(all_projs.shape[1]):
        a = array((all_projs[ii,jj,-3], all_projs[ii,jj,-2], all_projs[ii,jj,-1]))
        dist[ii,jj] = linalg.norm(a-b)



color_dimension = Z # change to desired fourth dimension
minn, maxx = color_dimension.min(), color_dimension.max()
norm = matplotlib.colors.LogNorm(minn, maxx)
m = plt.cm.ScalarMappable(norm=norm, cmap='Greys')
m.set_array([])
fcolors = m.to_rgba(color_dimension )

# Plot surface
ax3d.plot_surface(all_projs[:,:,-3], all_projs[:,:,-2], all_projs[:,:,-1], rstride=1, cstride=1, facecolors=fcolors,
                                                                            vmin=minn, vmax=maxx,alpha=1, shade=False,linewidth=0, antialiased=False)

# Plot best fit
ax3d.plot([popt_proj[-3]], [popt_proj[-2]], [popt_proj[-1]], 'o', c='#fb8072', ms=3, mec='k', mew=0.1, zorder=20)

# Plot data point
#ax3d.plot([data_proj[-3]], [data_proj[-2]], [data_proj[-1]], '*', c='#4daf4a', ms=5, zorder=20)

# Plot geodesics

ax3d.plot(all_gprojs[0][:,-3], all_gprojs[0][:,-2], all_gprojs[0][:,-1], '-', c='#bebada', lw=1, zorder=10)
ax3d.plot(all_gprojs[1][:,-3], all_gprojs[1][:,-2], all_gprojs[1][:,-1], '-', c='#80b1d3', lw=1, zorder=10)

# Plot other two points
ax3d.plot([p20_proj[-3]], [p20_proj[-2]], [p20_proj[-1]], 'o', ms=3, c='#fdb462', mec='k', mew=0.1, zorder=10)
ax3d.plot([p21_proj[-3]], [p21_proj[-2]], [p21_proj[-1]], 'o', ms=3, c='#b3de69', mec='k', mew=0.1, zorder=10)
ax3d.plot([p22_proj[-3]], [p22_proj[-2]], [p22_proj[-1]], 'o', ms=3, c='#ffffb3', mec='k', mew=0.1, zorder=10)

#plt.colorbar(m)

# Plot edges in black
ax3d.plot(all_projs[:,-1,-3], all_projs[:,-1,-2],
          all_projs[:,-1,-1], '-k', linewidth=0.1)
ax3d.plot(all_projs[:,0,-3], all_projs[:,0,-2],
          all_projs[:,0,-1], '-k', linewidth=0.1)
ax3d.plot(all_projs[0,:,-3], all_projs[0,:,-2],
          all_projs[0,:,-1], '-k', linewidth=0.1)
ax3d.plot(all_projs[-1,:,-3], all_projs[-1,:,-2],
          all_projs[-1,:,-1], '-k', linewidth=0.1)





ax3d.view_init(60, 130)
ax3d.set_xticks([])
ax3d.set_yticks([])
ax3d.set_zticks([])
ax3d.set_axis_off()


fig.set_size_inches(width, height)
fig.savefig('manifold_PCA.svg')

plt.show()







