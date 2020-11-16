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
#dipStartT = float(sys.argv[2])
#dipStopT = float(sys.argv[3])
#integralStartT = float(sys.argv[4])
#integralStopT = float(sys.argv[5])

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
               ('eStart', np.float64, (1)),
               ('deathTime', np.float64, (1)),
               ('rLenBack', np.uint32, (1))])

data = np.fromfile(fname, dtype=dt)

cut = np.bitwise_and(data['time'] < data['deathTime'], data['time'] > 0)

dataCut = data[cut]

print(len(dataCut))

plt.clf()
plt.hist(dataCut['time'], bins=np.arange(0,1700,1))
plt.show()

plt.clf()
plt.hist(data['zoff'], bins=20)
plt.show()

#zetas = []
#for row in dataCut:
#    majR = 0.5 if row['x'] < 0 else 1.0
#    minr = 1.0 if row['x'] < 0 else 0.5
#    zetas.append(minr - math.sqrt(row['x']**2 + (math.sqrt(row['y']**2 + (row['z']-row['zoff'])**2) - majR)**2))
#
#plt.hist(zetas, bins=100)
#plt.show()

print("Read "+str(len(data))+" Evts")

print(len(dataCut), np.mean(dataCut['time']))
