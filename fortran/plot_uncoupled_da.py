import matplotlib.pyplot as plt
import numpy

labels = ['Atmosphere', 'Ocean']
ranges = [(1, 21), (21, 37)]

for i, component in enumerate(['atm', 'ocn']):
    lyap = open('mean_lyapunov_{component}_3D.dat'.format(component=component)).read()
    lyap = numpy.array(lyap.split()[ranges[i][0]:ranges[i][1]], dtype='float')
    lyap *= 86400*1.032e-4
    length = (ranges[i][1] - ranges[i][0])
    plt.scatter(numpy.array([i*2.5]*length) + (numpy.arange(-length/2, length/2) + 0.5)/15, lyap, c='black', s=15)

plt.xticks([0, 2.5], labels)
plt.axhline(y=0.02,color='r',linestyle=':')
plt.axhline(y=-0.02,color='r',linestyle=':')
plt.ylabel('Lyapunov exponent (days$^{-1}$)')
plt.savefig('uncoupled_3D.pdf')
plt.clf()
