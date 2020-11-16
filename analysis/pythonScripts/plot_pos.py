#!/usr/bin/python3

import sys
import numpy as np
import matplotlib.pyplot as plt

jtonev = 6.2415091e27
np.random.seed(2303616184)

if(len(sys.argv) != 2):
	sys.exit("Error! Usage: ./plotArrivalTime fname")

fname = sys.argv[1]

dt = np.dtype([('rLenFront', np.uint32, (1)),
               ('energy', np.float64, (1)),
               ('theta', np.float64, (1)),
               ('time', np.float64, (1)),
               ('eperp', np.float64, (1)),
               ('x', np.float64, (1)),
               ('y', np.float64, (1)),
               ('z', np.float64, (1)),
               ('zoff', np.float64, (1)),
               ('nhit', np.int32, (1)),
               ('nhitBot', np.int32, (1)),
               ('nhitTop', np.int32, (1)),
               ('rLenBack', np.uint32, (1))])

data = np.fromfile(fname, dtype=dt)

#dataCut = data[np.logical_and(data['time'] < 80, data['time'] > 40)]
dataCut = data[data['time'] < 2000]
    
zeta = np.where(dataCut['x'] > 0,
                0.5 - np.sqrt(dataCut['x']**2 + (np.fabs(dataCut['z'] - dataCut['zoff']) - 1.0)**2),
                1.0 - np.sqrt(dataCut['x']**2 + (np.fabs(dataCut['z'] - dataCut['zoff']) - 0.5)**2))

plt.clf()
plt.hist(dataCut['time'], bins=2000)
plt.show()

plt.clf()
plt.hist(zeta, bins=200)
plt.show()

print(np.mean(dataCut['nhit']))
plt.clf()
plt.hist(dataCut['nhit'], bins=np.arange(0,50,1))
plt.show()

plt.clf()
plt.hist2d(dataCut['time'], 100*zeta, bins=100)
plt.xlabel('time [s]')
plt.ylabel('$\zeta$ [cm]')
plt.title("3 Dip Zeta by Dip")
plt.show()
