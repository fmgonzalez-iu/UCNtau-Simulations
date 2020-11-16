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

plt.hist(data['theta'], bins=100)
plt.show()