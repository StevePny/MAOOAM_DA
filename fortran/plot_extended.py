import numpy
import matplotlib.pyplot as plt

lyap_all = open('mean_lyapunov_297.dat').read()
lyap2_all = numpy.array(lyap_all.split()[1:37], dtype='float')

lyap_ocn = open('mean_lyapunov_299.dat').read()
lyap2_ocn = numpy.array(lyap_ocn.split()[1:53], dtype='float')

lyap_atm = open('mean_lyapunov_298.dat').read()
lyap2_atm = numpy.array(lyap_atm.split()[1:57], dtype='float')

for le in lyap2_all:
    lyap2_atm[numpy.argmin(abs(lyap2_atm - le))] = numpy.inf
    lyap2_ocn[numpy.argmin(abs(lyap2_ocn - le))] = numpy.inf

plt.scatter(range(36), lyap2_all)
plt.scatter(range(38, 38+20), lyap2_atm[numpy.isfinite(lyap2_atm)])
plt.scatter(range(38+20+2, 38+20+2+16), lyap2_ocn[numpy.isfinite(lyap2_ocn)])
