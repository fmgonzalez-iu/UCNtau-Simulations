#!/usr/bin/python3

import sys
import numpy as np
import matplotlib.pyplot as plt

jtonev = 6.2415091e27

def sumCts(data, dipStartT, dipStopT):
    integral = 0
    
    cutDip = (data['energy'] > 5.875/jtonev) & (data['time'] < data['deathTime']) & (data['time'] >= dipStartT) & (data['time'] <= dipStopT)
    weightsDip = ((data[cutDip]['energy']*jtonev)**1.216666667)*(np.cos(data[cutDip]['theta'])**0.25)
    
    dipSum = np.sum(weightsDip)
    dipSumErr = np.sqrt(np.sum(weightsDip*weightsDip))

    return (dipSum, dipSumErr)

def arrival(data, dipStartT, dipStopT):
    integral = 0
    sumT = 0
    
    cutDip = (data['energy'] > 5.875/jtonev) & (data['time'] < data['deathTime']) & (data['time'] >= dipStartT) & (data['time'] <= dipStopT)
    weightsDip = ((data[cutDip]['energy']*jtonev)**1.216666667)*(np.cos(data[cutDip]['theta'])**0.25)
    weightsArrival = data[cutDip]['time']*((data[cutDip]['energy']*jtonev)**1.216666667)*(np.cos(data[cutDip]['theta'])**0.25)
    
    return np.sum(weightsArrival)/np.sum(weightsDip)

np.random.seed(2303616184)

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

holdTs = [20, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]

data50 = np.fromfile("../datafiles/pse/pse_512k_50.bin", dtype=dt)
data1000 = np.fromfile("../datafiles/pse/pse_512k_1000.bin", dtype=dt)

(sum50, err50) = sumCts(data50, 500 + 50 + 20 + 20, 500 + 50 + 20 + 120)
(sum1000, err1000) = sumCts(data1000, 500 + 1000 + 20 + 20, 500 + 1000 + 20 + 120)

arrival50 = arrival(data50, 500 + 50 + 20 + 20, 500 + 50 + 20 + 120)
arrival1000 = arrival(data1000, 500 + 1000 + 20 + 20, 500 + 1000 + 20 + 120)

print((arrival1000-arrival50)/(np.log(sum50/sum1000)))

