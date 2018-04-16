import matplotlib.pyplot as plt
import numpy

labels = ['Coupled', 'Atmosphere', 'Ocean']
ranges = [(1, 37), (1, 21), (21, 37)]

for i, component in enumerate(['both', 'atm', 'ocn']):
    lyap = open('uncoupled/mean_lyapunov_{component}_099.dat'.format(component=component)).read()
    lyap = numpy.array(lyap.split()[ranges[i][0]:ranges[i][1]], dtype='float')
    length = (ranges[i][1] - ranges[i][0])
    plt.scatter(numpy.array([i*1.5]*length) + (numpy.arange(-length/2, length/2) + 0.5)/30, lyap, c='black', s=5)

plt.ylim(-0.2, 0.2)
plt.xticks([0, 1.5, 3], labels)
plt.ylabel('Lyapunov exponent (days$^{-1}$)')
plt.savefig('uncoupled.png')
plt.savefig('uncoupled.eps')
