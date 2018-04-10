import matplotlib.pyplot as plt
import xarray
import numpy

dat = xarray.DataArray(numpy.reshape(range(100), (10, 10)),
                       coords=[numpy.linspace(0.001, 1, 10),
                               numpy.linspace(0.001, 1, 10)],
                       dims=['coupling_motion', 'coupling_thermo'])

# Full motion coupling, gradual thermo
sim_ids = dat[9, :].data
thermo = dat[9, :]['coupling_thermo'].data

for i, sim_id in enumerate(sim_ids):
    lyap = open('mean_lyapunov_{sim_id}.dat'.format(sim_id=str(sim_id).zfill(3))).read()
    lyap = numpy.array(lyap.split('\n')[0].split()[1:], dtype='float')
    plt.scatter(numpy.array([thermo[i]]*36) + numpy.random.rand(36)/30, lyap, c='black', s=5)

plt.ylim(-0.5, 0.5)
plt.xlabel('Thermodynamic coupling')
plt.ylabel('Lyapunov exponent (days$^{-1}$)')
plt.savefig('gradual_thermo.png')

# Full thermo coupling, gradual motion
sim_ids = dat[:, 9].data
thermo = dat[:, 9]['coupling_motion'].data

for i, sim_id in enumerate(sim_ids):
    lyap = open('mean_lyapunov_{sim_id}.dat'.format(sim_id=str(sim_id).zfill(3))).read()
    lyap = numpy.array(lyap.split('\n')[0].split()[1:], dtype='float')
    plt.scatter(numpy.array([thermo[i]]*36) + numpy.random.rand(36)/30, lyap, c='black', s=5)


plt.ylim(-0.5, 0.5)
plt.xlabel('Motion coupling')
plt.ylabel('Lyapunov exponent (days$^{-1}$)')
plt.savefig('gradual_motion.png')

mat = numpy.zeros((10, 10))
for i in range(10):
    for j in range(10):
        sim_id = dat[i, j].data
        lyap = open('mean_lyapunov_{sim_id}.dat'.format(sim_id=str(sim_id).zfill(3))).read()
        lyap = numpy.array(lyap.split('\n')[0].split()[1:], dtype='float')
        mat[i, j] = sum(lyap > 0)

plt.matshow(mat)
plt.xticks(range(10), numpy.linspace(0.001, 1, 10))
plt.yticks(range(10), numpy.linspace(0.001, 1, 10))
plt.colorbar()
plt.xlabel('Motion coupling')
plt.ylabel('Thermodynamic coupling')
plt.title('Number of positive Lyapunov exponents')
plt.savefig('gradual_both.png', bbox_inches='tight')
