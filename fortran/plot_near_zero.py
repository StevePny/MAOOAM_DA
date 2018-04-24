import numpy
import matplotlib.pyplot as plt
from arch.bootstrap import StationaryBootstrap

f = open('lyapunov_exponents_297.dat', 'rb')

field = numpy.fromfile(f, dtype='float64').reshape([-1, 36])

bs = StationaryBootstrap(100, field)

def func(y):
    return y.mean(axis=0)

ci = bs.conf_int(func, 1000)

lyap_ocn = open('uncoupled/mean_lyapunov_ocn_099.dat').read()
lyap2_ocn = numpy.array(lyap_ocn.split()[21:37], dtype='float')

plt.plot(range(18), field.mean(axis=0)[2:20])
plt.plot(range(1,17), lyap2_ocn)
plt.vlines(range(18), ci[0][2:20], ci[1][2:20])
plt.hlines([0], xmin=0, xmax=17)

plt.title('Near-zero Lyapunov exponents')
plt.ylabel('Time$^{-1}$ (95% bootstrap confidence interval)')
plt.tight_layout()
plt.savefig('near_zero.pdf')
