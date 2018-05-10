import numpy
import matplotlib.pyplot as plt
from arch.bootstrap import StationaryBootstrap
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri

rpy2.robjects.numpy2ri.activate()

r_source = robjects.r['source']
block_length = r_source('ppw.R')[0]

def mean(y):
    return y.mean(axis=0)

with open('lyapunov_exponents_both_297.dat', 'rb') as f:
    lyap_coupled = numpy.fromfile(f, dtype='float64').reshape([-1, 36])

with open('lyapunov_exponents_ocn_297.dat', 'rb') as f: 
    lyap_ocn = numpy.fromfile(f, dtype='float64').reshape([-1, 36])[:, 20:36]

lyap_coupled *= 86400*1.032e-4
lyap_ocn *= 86400*1.032e-4

block_length_coupled = numpy.array(block_length(lyap_coupled))[:, 0]

all_ci = []
for var in range(36):
    bootstrap = StationaryBootstrap(block_length_coupled[var],
                                    lyap_coupled[:, var])


    ci = bootstrap.conf_int(mean, 1000)
    all_ci.append(ci)

ci_coupled = numpy.hstack(all_ci)

block_length_ocn = numpy.array(block_length(lyap_ocn))[:, 0]

all_ci = []
for var in range(16):
    bootstrap = StationaryBootstrap(block_length_ocn[var],
                                    lyap_ocn[:, var])


    ci = bootstrap.conf_int(mean, 1000)
    all_ci.append(ci)

ci_ocn = numpy.hstack(all_ci)

plt.scatter(range(17), lyap_coupled.mean(axis=0)[2:19], color='orange',
            s=5, zorder=1)
plt.scatter(numpy.arange(16) + 0.2, lyap_ocn.mean(axis=0), color='blue', s=5,
            zorder=2)

plt.vlines(numpy.arange(16) + 0.2, ci_ocn[0], ci_ocn[1], colors='blue')
plt.vlines(range(17), ci_coupled[0][2:19], ci_coupled[1][2:19],
           colors='orange')

plt.title('Near-zero Lyapunov exponents')
plt.ylabel('Days$^{-1}$ (95% bootstrap confidence interval)')

plt.hlines([0], xmin=0, xmax=16, linestyles='dashed', linewidth=0.3,
           zorder=0)

plt.tick_params(axis='x', which='both', bottom=False, top=False,
                labelbottom=False)

plt.legend(['Near-zero coupled LEs', 'Uncoupled ocean LEs'])
plt.tight_layout()
plt.savefig('near_zero.pdf')
