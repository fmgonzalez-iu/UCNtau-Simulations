#!/usr/bin/python3

import sys
import numpy as np
import matplotlib.pyplot as plt

jtonev = 6.2415091e27

def sumCts(data, dipStartT, dipStopT, integralStartT, integralStopT):
    integral = 0
    norm = 0
    cutDip = (data['energy'] > 5.875/jtonev) & (data['time'] < data['deathTime'] + 150) & (data['time'] >= dipStartT) & (data['time'] <= dipStopT)
    cutNorm = (data['energy'] > 5.875/jtonev) & (data['time'] < data['deathTime'] + 150) & (data['time'] >= integralStartT) & (data['time'] <= integralStopT)
    weightsDip = ((data[cutDip]['energy']*jtonev)**1.216666667)*(np.cos(data[cutDip]['theta'])**0.25)
    weightsNorm = ((data[cutNorm]['energy']*jtonev)**1.216666667)*(np.cos(data[cutNorm]['theta'])**0.25)
    

    dipSum = np.sum(weightsDip)
    norm = np.sum(weightsNorm)
    dipSumErr = np.sqrt(np.sum(weightsDip*weightsDip))

    return (dipSum/norm, dipSumErr/norm)

def sumPSE(data, holdT):
    tStartDip = 350 + holdT + 20
    tEndDip = 350 + holdT + 20 + 20
    tEndInt = 350 + holdT + 20 + 100
    return sumCts(data, tStartDip, tEndDip, tStartDip, tEndInt)

def plotPSE(data, holdT):
    tStartDip = 350 + holdT + 20
    tEndDip = 350 + holdT + 20 + 20
    tEndInt = 350 + holdT + 20 + 100
    cutNorm = (data['energy'] > 5.875/jtonev) & (data['time'] < data['deathTime']) & (data['time'] >= tStartDip) & (data['time'] <= tEndInt)
    weightsNorm = ((data[cutNorm]['energy']*jtonev)**1.216666667)*(np.cos(data[cutNorm]['theta'])**0.25)
    plt.clf()
    plt.hist(data[cutNorm]['time'], bins=np.arange(tStartDip,tEndInt+1,1), weights=weightsNorm)
    plt.hist(data[cutNorm]['time'], bins=np.arange(tStartDip,tEndDip+1,1), weights=weightsNorm)
    plt.show()
    
def plotPos(data, holdT):
    tStartDip = 350 + holdT + 20
    tEndDip = 350 + holdT + 20 + 20
    tEndInt = 350 + holdT + 20 + 100
    cutNorm = (data['energy'] > 5.875/jtonev) & (data['time'] < data['deathTime']) & (data['time'] >= tStartDip) & (data['time'] <= tEndDip) & (data['energy'] > 35/jtonev)
    weightsNorm = ((data[cutNorm]['energy']*jtonev)**1.216666667)*(np.cos(data[cutNorm]['theta'])**0.25)
    zetas = (1.0 - 0.5/(1 + np.exp(-1000*data[cutNorm]['x']))
             - np.sqrt(
                 data[cutNorm]['x']**2
                 +
                 (
                     np.sqrt(data[cutNorm]['y']**2 + (data[cutNorm]['z']-data[cutNorm]['zoff'])**2)
                     - 0.5 -0.5/(1 + np.exp(-1000*data[cutNorm]['x'])))
                 **2)
            )
    print(np.mean(zetas))
    plt.clf()
    plt.hist(zetas, bins=100, weights=weightsNorm)
    plt.show()
    
def plotSpect(data, holdT):
    tStartDip = 350 + holdT + 20
    tEndDip = 350 + holdT + 20 + 20
    tEndInt = 350 + holdT + 20 + 100
    cutNorm = (data['energy'] > 5.875/jtonev) & (data['time'] < data['deathTime']) & (data['time'] >= tStartDip) & (data['time'] <= tEndInt)
    weightsNorm = ((data[cutNorm]['energy']*jtonev)**1.216666667)*(np.cos(data[cutNorm]['theta'])**0.25)
    plt.clf()
    plt.hist(data[cutNorm]['energy']*jtonev, bins=100, weights=weightsNorm)
    plt.show()
    plt.clf()
    plt.hist(data[cutNorm]['theta'], bins=100, weights=weightsNorm)
    plt.show()
    plt.clf()
    plt.hist(data['energy']*jtonev, bins=100)
    plt.show()
    plt.clf()
    plt.hist(data['theta'], bins=100)
    plt.show()

def plotEDiff(data, holdT):
    tStartDip = 350 + holdT + 20
    tEndDip = 350 + holdT + 20 + 20
    tEndInt = 350 + holdT + 20 + 100
    cutNorm = (data['energy'] > 5.875/jtonev) & (data['time'] < data['deathTime']) & (data['time'] >= tStartDip) & (data['time'] <= tEndInt)
    weightsNorm = ((data[cutNorm]['energy']*jtonev)**1.216666667)*(np.cos(data[cutNorm]['theta'])**0.25)
    plt.clf()
    plt.hist((data[cutNorm]['energy'] - data[cutNorm]['eStart'])*jtonev, bins=100, weights=weightsNorm)
    plt.yscale('log')
    plt.show()
    print(np.sum(weightsNorm*(data[cutNorm]['energy'] - data[cutNorm]['eStart'])*jtonev)/np.sum(weightsNorm))

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
#holdTs = [20, 500, 1000]
fNames = {}
for t in holdTs:
    fNames[t] = "../datafiles/pse/pse_heat_512k_"+str(t)+".bin"
#    fNames[t] = "../datafiles/pse/pse_6ms_512k_"+str(t)+".bin"

#data = np.fromfile(fNames[20], dtype=dt)
#(shortCts, shortErr) = sumCts(data, 200+20, 200+20+120, 200+20+20, 200+20+40)

for t in holdTs:
    data = np.fromfile(fNames[t], dtype=dt)
    assert(np.all(data['rLenFront']==data[0]['rLenFront']))
#    plotPSE(data, t)
#    plotPos(data, t)
#    plotSpect(data, t)
#    plotEDiff(data, t)
    (frac, fracErr) = sumPSE(data, t)
    print(t, 100*frac, 0.0, 100*fracErr)
    
    #(cts, err) = sumCts(data, 200+t, 200+t+120, 200+t+20, 200+t+40)
    #print(t, cts/shortCts, 0.0, err/shortCts)
