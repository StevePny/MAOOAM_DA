import matplotlib.pyplot as plt
import xarray
import numpy

dat = xarray.open_dataarray('sim_ids.nc')

# Full motion coupling, gradual thermo
sim_ids = dat[9, :].data
thermo = dat[9, :]['coupling_thermo'].data

for i, sim_id in enumerate(sim_ids):
    lyap = open('mean_lyapunov_{sim_id}.dat'.format(sim_id=sim_id)).read()
    lyap = numpy.array(lyap.split('\n')[0].split()[1:], dtype='float')
    plt.scatter(numpy.array([thermo[i]]*36) + numpy.random.rand(36)/30, lyap, c='black', s=5)

plt.ylim(-0.2, 0.2)
plt.xlabel('Thermodynamic coupling')
plt.ylabel('Lyapunov exponent (days$^{-1}$)')
plt.savefig('gradual_thermo.png')
plt.clf()

# Full thermo coupling, gradual motion
sim_ids = dat[:, 9].data
thermo = dat[:, 9]['coupling_motion'].data

for i, sim_id in enumerate(sim_ids):
    lyap = open('mean_lyapunov_{sim_id}.dat'.format(sim_id=sim_id)).read()
    lyap = numpy.array(lyap.split('\n')[0].split()[1:], dtype='float')
    plt.scatter(numpy.array([thermo[i]]*36) + numpy.random.rand(36)/30, lyap, c='black', s=5)

plt.ylim(-0.2, 0.2)
plt.xlabel('Motion coupling')
plt.ylabel('Lyapunov exponent (days$^{-1}$)')
plt.savefig('gradual_motion.png')
plt.clf()

# Gradual thermo, gradual motion
sim_ids = numpy.diag(dat)
thermo = dat[:, 9]['coupling_motion'].data

for i, sim_id in enumerate(sim_ids):
    lyap = open('mean_lyapunov_{sim_id}.dat'.format(sim_id=sim_id)).read()
    lyap = numpy.array(lyap.split('\n')[0].split()[1:], dtype='float')
    plt.scatter(numpy.array([thermo[i]]*36) + numpy.random.rand(36)/30, lyap, c='black', s=5)

plt.ylim(-0.2, 0.2)
plt.xlabel('Coupling')
plt.ylabel('Lyapunov exponent (days$^{-1}$)')
plt.savefig('gradual_both.png')
plt.clf()

# Low thermo, gradual motion
sim_ids = dat[:, 0].data
thermo = dat[:, 0]['coupling_motion'].data

for i, sim_id in enumerate(sim_ids):
    lyap = open('mean_lyapunov_{sim_id}.dat'.format(sim_id=sim_id)).read()
    lyap = numpy.array(lyap.split('\n')[0].split()[1:], dtype='float')
    plt.scatter(numpy.array([thermo[i]]*36) + numpy.random.rand(36)/30, lyap, c='black', s=5)

plt.ylim(-0.2, 0.2)
plt.xlabel('Motion coupling')
plt.ylabel('Lyapunov exponent (days$^{-1}$)')
plt.savefig('gradual_low_thermo.png')
plt.clf()

# Low motion, gradual thermo
sim_ids = dat[0, :].data
thermo = dat[0, :]['coupling_thermo'].data

for i, sim_id in enumerate(sim_ids):
    lyap = open('mean_lyapunov_{sim_id}.dat'.format(sim_id=sim_id)).read()
    lyap = numpy.array(lyap.split('\n')[0].split()[1:], dtype='float')
    plt.scatter(numpy.array([thermo[i]]*36) + numpy.random.rand(36)/30, lyap, c='black', s=5)

plt.ylim(-0.2, 0.2)
plt.xlabel('Thermodynamic coupling')
plt.ylabel('Lyapunov exponent (days$^{-1}$)')
plt.savefig('gradual_low_motion.png')
plt.clf()

# Plot number of positive LEs
mat = numpy.zeros((10, 10), dtype='int')
for i in range(10):
    for j in range(10):
        sim_id = dat[i, j].data
        lyap = open('mean_lyapunov_{sim_id}.dat'.format(sim_id=sim_id)).read()
        lyap = numpy.array(lyap.split('\n')[0].split()[1:], dtype='float')
        mat[i, j] = sum(lyap > 0)

N = mat.max() - mat.min() + 1
plt.matshow(mat, fignum=1, origin='lower', cmap=plt.cm.get_cmap('viridis', N))
ax = plt.gca()
ax.xaxis.tick_bottom()
plt.xticks(range(10), numpy.around(numpy.linspace(0.001, 1, 10), decimals=2))
plt.yticks(range(10), numpy.around(numpy.linspace(0.001, 1, 10), decimals=2))
plt.colorbar(ticks=range(mat.min(), mat.max() + 1))
plt.clim(mat.min() - 0.5, mat.max() + 0.5)
plt.xlabel('Motion coupling')
plt.ylabel('Thermodynamic coupling')
plt.title('Number of positive Lyapunov exponents')
plt.savefig('gradual_grid.png', bbox_inches='tight')
