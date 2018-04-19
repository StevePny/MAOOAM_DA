import matplotlib.pyplot as plt
import xarray
import numpy

dat = xarray.open_dataarray('sim_ids.nc')

# Full momentum coupling, gradual thermo
sim_ids = dat[9, :].data
thermo = dat[9, :]['coupling_thermo'].data

for i, sim_id in enumerate(sim_ids):
    lyap = open('coupled/mean_lyapunov_{sim_id}.dat'.format(sim_id=sim_id)).read()
    lyap = numpy.array(lyap.split('\n')[0].split()[1:], dtype='float')
    lyap *= 86400*1.032e-4
    plt.scatter(numpy.array([thermo[i]]*36) + numpy.random.rand(36)/30, lyap, c='black', s=5)

plt.ylim(-1.8, 1.8)
plt.xlabel('Thermodynamic coupling')
plt.ylabel('Lyapunov exponent (days$^{-1}$)')
plt.savefig('gradual_thermo.pdf')
plt.clf()

# Full thermo coupling, gradual momentum
sim_ids = dat[:, 9].data
thermo = dat[:, 9]['coupling_motion'].data

for i, sim_id in enumerate(sim_ids):
    lyap = open('coupled/mean_lyapunov_{sim_id}.dat'.format(sim_id=sim_id)).read()
    lyap = numpy.array(lyap.split('\n')[0].split()[1:], dtype='float')
    lyap *= 86400*1.032e-4
    plt.scatter(numpy.array([thermo[i]]*36) + numpy.random.rand(36)/30, lyap, c='black', s=5)

plt.ylim(-1.8, 1.8)
plt.xlabel('Momentum coupling')
plt.ylabel('Lyapunov exponent (days$^{-1}$)')
plt.savefig('gradual_momentum.pdf')
plt.clf()

# Gradual thermo, gradual momentum
sim_ids = numpy.diag(dat)
thermo = dat[:, 9]['coupling_motion'].data

for i, sim_id in enumerate(sim_ids):
    lyap = open('coupled/mean_lyapunov_{sim_id}.dat'.format(sim_id=sim_id)).read()
    lyap = numpy.array(lyap.split('\n')[0].split()[1:], dtype='float')
    lyap *= 86400*1.032e-4
    plt.scatter(numpy.array([thermo[i]]*36) + numpy.random.rand(36)/30, lyap, c='black', s=5)

plt.ylim(-1.8, 1.8)
plt.xlabel('Coupling')
plt.ylabel('Lyapunov exponent (days$^{-1}$)')
plt.savefig('gradual_both.pdf')
plt.clf()

# Low thermo, gradual momentum
sim_ids = dat[:, 0].data
thermo = dat[:, 0]['coupling_motion'].data

for i, sim_id in enumerate(sim_ids):
    lyap = open('coupled/mean_lyapunov_{sim_id}.dat'.format(sim_id=sim_id)).read()
    lyap = numpy.array(lyap.split('\n')[0].split()[1:], dtype='float')
    lyap *= 86400*1.032e-4
    plt.scatter(numpy.array([thermo[i]]*36) + numpy.random.rand(36)/30, lyap, c='black', s=5)

plt.ylim(-1.8, 1.8)
plt.xlabel('Momentum coupling')
plt.ylabel('Lyapunov exponent (days$^{-1}$)')
plt.savefig('gradual_low_thermo.pdf')
plt.clf()

# Low momentum, gradual thermo
sim_ids = dat[0, :].data
thermo = dat[0, :]['coupling_thermo'].data

for i, sim_id in enumerate(sim_ids):
    lyap = open('coupled/mean_lyapunov_{sim_id}.dat'.format(sim_id=sim_id)).read()
    lyap = numpy.array(lyap.split('\n')[0].split()[1:], dtype='float')
    lyap *= 86400*1.032e-4
    plt.scatter(numpy.array([thermo[i]]*36) + numpy.random.rand(36)/30, lyap, c='black', s=5)

plt.ylim(-1.8, 1.8)
plt.xlabel('Thermodynamic coupling')
plt.ylabel('Lyapunov exponent (days$^{-1}$)')
plt.savefig('gradual_low_momentum.pdf')
plt.clf()

# Plot number of positive LEs
mat = numpy.zeros((10, 10), dtype='int')
for i in range(10):
    for j in range(10):
        sim_id = dat[i, j].data
        lyap = open('coupled/mean_lyapunov_{sim_id}.dat'.format(sim_id=sim_id)).read()
        lyap = numpy.array(lyap.split('\n')[0].split()[1:], dtype='float')
        lyap *= 86400*1.032e-4
        mat[i, j] = sum(lyap > 0)

N = mat.max() - mat.min() + 1
plt.matshow(mat, fignum=1, origin='lower', cmap=plt.cm.get_cmap('viridis', N))
ax = plt.gca()
ax.xaxis.tick_bottom()
plt.xticks(range(10), numpy.around(numpy.linspace(0.001, 1, 10), decimals=2))
plt.yticks(range(10), numpy.around(numpy.linspace(0.001, 1, 10), decimals=2))
plt.colorbar(ticks=range(mat.min(), mat.max() + 1))
plt.clim(mat.min() - 0.5, mat.max() + 0.5)
plt.xlabel('Momentum coupling')
plt.ylabel('Thermodynamic coupling')
plt.title('Number of positive Lyapunov exponents')
plt.savefig('gradual_grid.pdf', bbox_inches='tight')
plt.clf()

# Plot number of positive and neutral LEs
mat = numpy.zeros((10, 10), dtype='int')
for i in range(10):
    for j in range(10):
        sim_id = dat[i, j].data
        lyap = open('coupled/mean_lyapunov_{sim_id}.dat'.format(sim_id=sim_id)).read()
        lyap = numpy.array(lyap.split('\n')[0].split()[1:], dtype='float')
        lyap *= 86400*1.032e-4
        mat[i, j] = sum(lyap > -0.02)

N = mat.max() - mat.min() + 1
plt.matshow(mat, fignum=1, origin='lower', cmap=plt.cm.get_cmap('viridis', N))
ax = plt.gca()
ax.xaxis.tick_bottom()
plt.xticks(range(10), numpy.around(numpy.linspace(0.001, 1, 10), decimals=2))
plt.yticks(range(10), numpy.around(numpy.linspace(0.001, 1, 10), decimals=2))
plt.colorbar(ticks=range(mat.min(), mat.max() + 1))
plt.clim(mat.min() - 0.5, mat.max() + 0.5)
plt.xlabel('Momentum coupling')
plt.ylabel('Thermodynamic coupling')
plt.title('Number of positive and neutral Lyapunov exponents')
plt.savefig('gradual_grid_neutral.pdf', bbox_inches='tight')
