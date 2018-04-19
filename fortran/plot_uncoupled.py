import matplotlib.pyplot as plt
import numpy

labels = ['Coupled', 'Atmosphere', 'Ocean']
ranges = [(1, 37), (1, 21), (21, 37)]

plt.figure(figsize=(12, 5))

for i, component in enumerate(['both', 'atm', 'ocn']):
    lyap = open('uncoupled/mean_lyapunov_{component}_099.dat'.format(component=component)).read()
    lyap = numpy.array(lyap.split()[ranges[i][0]:ranges[i][1]], dtype='float')
    lyap *= 86400*1.032e-4
    length = (ranges[i][1] - ranges[i][0])
    plt.scatter(numpy.array([i*2.5]*length) + (numpy.arange(-length/2, length/2) + 0.5)/15, lyap, c='black', s=15)

plt.ylim(-1, 0.5)
plt.xticks([0, 2.5, 5], labels)
plt.axhline(y=0.02,color='r',linestyle=':')
plt.axhline(y=-0.02,color='r',linestyle=':')
plt.ylabel('Lyapunov exponent (days$^{-1}$)')
plt.savefig('uncoupled.png')
plt.savefig('uncoupled.pdf')
plt.clf()

plt.figure(figsize=(12, 5))

for i, component in enumerate(['both', 'atm', 'ocn']):
    lyap = open('uncoupled/mean_lyapunov_{component}_stochastic_099.dat'.format(component=component)).read()
    lyap = numpy.array(lyap.split()[ranges[i][0]:ranges[i][1]], dtype='float')
    lyap *= 86400*1.032e-4
    length = (ranges[i][1] - ranges[i][0])
    plt.scatter(numpy.array([i*2.5]*length) + (numpy.arange(-length/2, length/2) + 0.5)/15, lyap, c='black', s=15)

plt.ylim(-1, 0.5)
plt.xticks([0, 2.5, 5], labels)
plt.axhline(y=0.02,color='r',linestyle=':')
plt.axhline(y=-0.02,color='r',linestyle=':')
plt.ylabel('Lyapunov exponent (days$^{-1}$)')
plt.savefig('uncoupled_stochastic.png')
plt.savefig('uncoupled_stochastic.pdf')
