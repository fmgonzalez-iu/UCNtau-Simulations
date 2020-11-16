#!/usr/bin/python3

import sys
import numpy as np
import math
import matplotlib.pyplot as plt

jtonev = 6.2415091e27
np.random.seed(2303616184)

if(len(sys.argv) != 2):
	sys.exit("Error! Usage: ./plotArrivalTime fname")

fname = sys.argv[1]

dt = np.dtype([('rLenFront', np.uint32, (1)),
               ('energy', np.float64, (1)),
               ('theta', np.float64, (1)),
               ('times', np.float32, (50)),
               ('eperps', np.float32, (50)),
               ('rLenBack', np.uint32, (1))])

data = np.fromfile(fname, dtype=dt)

#wgts = np.ones(np.shape(data[data['energy'] > 7/jtonev]['eperps']))
##wgts = wgts * jtonev * data[data['energy'] > 7/jtonev]['energy'][:,np.newaxis]
#plt.clf()
#plt.hist(data[data['energy'] > 7/jtonev]['eperps'].flatten()*jtonev, bins=np.arange(1,40,1), weights=wgts.flatten())
#plt.show()

#wgts = np.ones(np.shape(data[data['energy'] > 7/jtonev]['times']))
#wgts = wgts * jtonev * data[data['energy'] > 7/jtonev]['energy'][:,np.newaxis]
#plt.clf()
##plt.hist(data[data['energy'] > 7/jtonev]['times'].flatten(), bins=np.arange(300+20+41, 300+20+41+184, 1), weights=wgts.flatten())
#plt.hist(data[data['energy'] > 7/jtonev]['times'].flatten(), bins=np.arange(41, 41+184, 1), weights=wgts.flatten())
#plt.show()

#wgts = np.ones(np.shape(data[data['energy'] > 7/jtonev]['energy']))
#wgts = wgts * (jtonev * data[data['energy'] > 7/jtonev]['energy'])**(1.227273-1)
#plt.clf()
#plt.hist(data[data['energy'] > 7/jtonev]['energy']*jtonev, bins=np.arange(1,40,1), weights=wgts)
#plt.show()

#plt.clf()
#plt.hist(data['theta'], bins=100)
#plt.show()

#xs = []
#ys = []
#for e in np.arange(7,30,0.1):
#    cut1 = np.bitwise_and(data['energy'] > e/jtonev, data['energy'] < (e+0.1)/jtonev)
#    cut2 = np.bitwise_and(data['theta'] > np.pi/4-np.pi/100, data['theta'] < np.pi/4+np.pi/100)
#    cut = np.bitwise_and(cut1, cut2)
#    xs.append(e)
##    ys.append(np.mean(data[cut]['eperps'].flatten()))
#    ys.append(np.sum(data[cut]['eperps'].flatten())/(49*len(data[cut])))
#plt.plot(xs, ys)
#plt.show()

#xs = []
#ys = []
#cut1 = np.bitwise_and(data['energy'] > 20/jtonev, data['energy'] < (20+0.1)/jtonev)
#cut2 = np.bitwise_and(data['theta'] > np.pi/4-np.pi/100, data['theta'] < np.pi/4+np.pi/100)
#cut = np.bitwise_and(cut1, cut2)
#for i in range(0, 50):
#    print(i, np.mean(data[cut]['eperps'][:,i]))

plt.clf()
plt.hist2d(data['energy'], data['theta'], bins=100)
plt.show()