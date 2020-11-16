#!/usr/bin/python3

import sys
import numpy as np
import math
#import matplotlib.pyplot as plt

jtonev = 6.2415091e27
np.random.seed(2303616184)

if(len(sys.argv) != 6):
	sys.exit("Error! Usage: ./plotArrivalTime fname dipStartT dipStopT integralStartT integralStopT")

fname = sys.argv[1]
dipStartT = float(sys.argv[2])
dipStopT = float(sys.argv[3])
integralStartT = float(sys.argv[4])
integralStopT = float(sys.argv[5])

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

print("Read "+str(len(data))+" Evts")

#plt.clf()
#plt.hist(data['theta'], bins=100)
#plt.show()
#plt.clf()
#plt.hist(data['energy']*jtonev, bins=100)
#plt.show()

integral = 0
intWgtSqr = 0
norm = 0
normWgtSqr = 0

def weight(row):
#    return 1
    if row['energy']*jtonev < 11.5:
        return 0
#    return math.pow(row['energy']*jtonev/34.5, 1.5)*math.pow(np.sin(row['theta']), 1)*math.pow(np.cos(row['theta']),(2))
    return math.pow(row['energy']*jtonev/34.5, 1.227273)*np.sin(row['theta'])*math.pow(np.cos(row['theta']),(1+0.181818)) #minimum
#    return math.pow(row['energy']*jtonev/34.5, 1.4)*np.sin(row['theta'])*math.pow(np.cos(row['theta']),(1+0.1))
#    return math.pow(row['energy']*jtonev/34.5, 1.75)*np.sin(row['theta'])*math.pow(np.cos(row['theta']),(1+0.5))
#    return math.pow(row['energy']*jtonev/34.5, 1.2)*np.sin(row['theta'])*math.pow(np.cos(row['theta']),(1+0.15))

meanTime = 0

for row in data:
    wgt = weight(row)
    if row['time'] >= dipStartT and row['time'] < dipStopT:
        norm += wgt
        normWgtSqr += wgt*wgt
        meanTime += wgt*row['time']
    if row['time'] >= integralStartT and row['time'] < integralStopT:
        integral += wgt
        intWgtSqr += wgt*wgt

meanTime = meanTime / norm
print(integral, np.sqrt(intWgtSqr), norm, np.sqrt(normWgtSqr), meanTime)


#print(integral/norm, np.sqrt(integral)/norm)
print(integral/norm, integral/norm*np.sqrt(intWgtSqr/(integral**2) + normWgtSqr/(norm**2)))
