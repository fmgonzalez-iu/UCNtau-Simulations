#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

jtonev = 6.2415091e27
mass_n = 1.674927471e-27
grav = 9.80665
minz = -1.464413669130002
minu = -2.390352484438862e-26
eCleanLow = 5.571749397933261e-27
eCleanHigh = 6.396044292436557e-27
    
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
    tStartDip = holdT + 20
    tEndDip = holdT + 20 + 20
    tEndInt = holdT + 20 + 120
    cutNorm = (data['time'] < data['deathTime']) & (data['time'] >= tStartDip) & (data['time'] <= tEndInt)
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
    tStartDip = holdT + 20
    tEndDip = holdT + 20 + 20
    tEndInt = holdT + 20 + 120
    cutNorm = (data['time'] < data['deathTime']) & (data['time'] >= tStartDip) & (data['time'] <= tEndInt)
    plt.clf()
    plt.hist((data[cutNorm]['energy'] - data[cutNorm]['eStart'])*jtonev, bins=100)
    plt.yscale('log')
    plt.show()

def countUCN(data, holdT):
    tStart = holdT + 20
    tEnd = holdT + 20 + 120
    cut = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEnd)
    return np.sum(cut)

def avgT(data, holdT):
    tStart = holdT + 20
    tEnd = holdT + 20 + 120
    cut = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEnd)
    return np.mean(data[cut]['time'])

def plotDip2(data, holdT):
    tStart = holdT + 20
    tEndDip2 = holdT + 20 + 20
    tEnd = holdT + 20 + 120
    cut = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEnd)
    cutDip2 = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEndDip2)
    plt.clf()
    plt.hist(data[cut]['time'], bins=np.arange(tStart, tEnd+1, 1))
    plt.hist(data[cutDip2]['time'], bins=np.arange(tStart, tEnd+1, 1), color='red')
    plt.show()

def plotDip1(data, holdT):
    tStart = 0.0
    tEndDip1 = holdT + 20
    tStartDip2 = holdT + 20
    tEndDip2 = holdT + 20 + 20
    tEnd = holdT + 20 + 120
    cut = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEnd)
    cutDip1 = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEndDip1)
    cutDip2 = (data['zoff'] > 0) & (data['time'] > tStartDip2) & (data['time'] < tEndDip2)
    plt.clf()
    plt.hist(data[cut]['time'], bins=np.arange(tStart, tEnd+1, 1), color='blue')
    plt.hist(data[cutDip2]['time'], bins=np.arange(tStart, tEnd+1, 1), color='red')
    plt.hist(data[cutDip1]['time'], bins=np.arange(tStart, tEnd+1, 1), color='green')
    plt.yscale('log')
    plt.show()

def plotDip1HighE(data, holdT):
    tStart = holdT
    tEndDip1 = holdT + 20
    tStartDip2 = holdT + 20
    tEndDip2 = holdT + 20 + 20
    tEnd = holdT + 20 + 120
    cut = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEnd) & (data['eStart'] >= eCleanLow) & (data['energy'] >= eCleanLow)
    cutAll = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEnd)
    cutDip1 = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEndDip1) & (data['eStart'] >= eCleanLow) & (data['energy'] >= eCleanLow)
#    cutDip2 = (data['zoff'] > 0) & (data['time'] > tStartDip2) & (data['time'] < tEndDip2) & (data['eStart'] >= eCleanLow) & (data['energy'] >= eCleanLow)
    plt.clf()
    plt.hist(data[cutAll]['time'], bins=np.arange(tStart, tEnd+1, 1), color='blue', label='All UCN')
    plt.hist(data[cut]['time'], bins=np.arange(tStart, tEnd+1, 1), color='red', label='E_i > 38cm & E_f > 38cm')
#    plt.hist(data[cutDip2]['time'], bins=np.arange(tStart, tEnd+1, 1), color='red')
#    plt.hist(data[cutDip1]['time'], bins=np.arange(tStart, tEnd+1, 1), color='green')
    plt.xlabel("Time [s]")
    plt.ylabel("Counts")
    plt.text(20, 50, r'$\Sigma$ = '+str(np.sum(cutDip1)))
    plt.text(80, 10, r'$\Sigma$ = '+str(np.sum(cut)))
    plt.title("Uncleaned Counting Efficiency")
    plt.yscale('log')
    plt.legend()
    plt.show()

def plotDip1Heat(data, holdT):
    tStart = holdT
    tEndDip1 = holdT + 20
    tStartDip2 = holdT + 20
    tEndDip2 = holdT + 20 + 20
    tEnd = holdT + 20 + 120
    cut = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEnd) & (data['eStart'] < eCleanLow) & (data['energy'] > eCleanLow)
    cutAll = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEnd)
    cutDip1 = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEndDip1) & (data['eStart'] < eCleanLow) & (data['energy'] > eCleanLow)
#    cutDip2 = (data['zoff'] > 0) & (data['time'] > tStartDip2) & (data['time'] < tEndDip2) & (data['eStart'] < eCleanLow) & (data['energy'] > eCleanLow)
    plt.clf()
    plt.hist(data[cutAll]['time'], bins=np.arange(tStart, tEnd+1, 1), color='blue', label='All UCN')
    plt.hist(data[cut]['time'], bins=np.arange(tStart, tEnd+1, 1), color='red', label='E_i < 38cm & E_f > 38cm')
#    plt.hist(data[cutDip2]['time'], bins=np.arange(tStart, tEnd+1, 1), color='red')
#    plt.hist(data[cutDip1]['time'], bins=np.arange(tStart, tEnd+1, 1), color='green')
    plt.xlabel("Time [s]")
    plt.ylabel("Counts")
    plt.text(1400, 50, r'$\Sigma$ = '+str(np.sum(cutDip1)))
    plt.text(1460, 10, r'$\Sigma$ = '+str(np.sum(cut)))
    plt.title("Heating Counting Efficiency")
    plt.yscale('log')
    plt.legend()
    plt.show()

def plotDip1_43cm(data, holdT):
    tStart = holdT
    tEndDip1 = holdT + 20
    tStartDip2 = holdT + 20
    tEndDip2 = holdT + 20 + 20
    tEnd = holdT + 20 + 120
    cut = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEnd) & (data['energy'] > eCleanHigh)
    cutAll = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEnd)
    cutDip1 = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEndDip1) & (data['energy'] > eCleanHigh)
#    cutDip2 = (data['zoff'] > 0) & (data['time'] > tStartDip2) & (data['time'] < tEndDip2) & (data['eStart'] < eCleanLow) & (data['energy'] > eCleanLow)
    plt.clf()
    plt.hist(data[cutAll]['time'], bins=np.arange(tStart, tEnd+1, 1), color='blue', label='All UCN')
    plt.hist(data[cut]['time'], bins=np.arange(tStart, tEnd+1, 1), color='red', label='E_f > 43cm')
#    plt.hist(data[cutDip2]['time'], bins=np.arange(tStart, tEnd+1, 1), color='red')
#    plt.hist(data[cutDip1]['time'], bins=np.arange(tStart, tEnd+1, 1), color='green')
    plt.xlabel("Time [s]")
    plt.ylabel("Counts")
    plt.text(1400, 50, r'$\Sigma$ = '+str(np.sum(cutDip1)))
    plt.text(1460, 10, r'$\Sigma$ = '+str(np.sum(cut)))
    plt.title("Heating Counting Efficiency")
    plt.yscale('log')
    plt.legend()
    plt.show()

def plotDip1LowE(data, holdT):
    tStart = holdT
    tEndDip1 = holdT + 20
    tStartDip2 = holdT + 20
    tEndDip2 = holdT + 20 + 20
    tEnd = holdT + 20 + 120
    cut = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEnd) & (data['eStart'] <= 0.2*mass_n*grav)
    cutDip1 = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEndDip1) & (data['eStart'] <= 0.2*mass_n*grav)
    cutDip2 = (data['zoff'] > 0) & (data['time'] > tStartDip2) & (data['time'] < tEndDip2) & (data['eStart'] <= 0.2*mass_n*grav)
    plt.clf()
    plt.hist(data[cut]['time'], bins=np.arange(tStart, tEnd+1, 1), color='blue')
    plt.hist(data[cutDip2]['time'], bins=np.arange(tStart, tEnd+1, 1), color='red')
    plt.hist(data[cutDip1]['time'], bins=np.arange(tStart, tEnd+1, 1), color='green')
    plt.yscale('log')
    plt.show()
    
def ratDip2(data, holdT):
    tStart = holdT + 20
    tEndDip2 = holdT + 20 + 20
    tEnd = holdT + 20 + 120
    cut = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEnd)
    cutDip2 = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEndDip2)
    sumDip1 = np.sum(cutDip2)
    norm = np.sum(cut)
    return (sumDip1/norm, sumDip1/norm*np.sqrt(1/sumDip1 + 1/norm))

def sumCleanChk(data, holdT):
    tStart = holdT
    tEnd = holdT + 20
    cut = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEnd)
    return np.sum(cut)

def plotClean(data, holdT):
    cut = (data['zoff'] == -2) & (data['time'] < 0)
    plt.clf()
#    plt.hist(data[cut]['time'] + data[cut]['settlingTime'], bins=np.arange(0,50,1))
#    plt.hist(data[cut]['time'] + 50, bins=np.arange(0,100,1))
#    plt.yscale('log')
#    plt.hist(data[cut]['energy']/(mass_n*grav)+(minz+1.5), bins=100)
#    plt.hist2d(data[cut]['energy']/(mass_n*grav)+(minz+1.5), bins=100)
#    plt.show()

def plotStorageClean(data, holdT):
    cutUncl = (data['zoff'] == -3) & (data['eStart'] >= eCleanLow)
    cutHeat = (data['zoff'] == -3) & (data['eStart'] < eCleanLow)
#    cutAll = (data['zoff'] == -3)
    plt.clf()
    plt.hist(data[cutHeat]['time'], bins=np.arange(0,holdT,1), color='blue', label='Heated')
    plt.hist(data[cutUncl]['time'], bins=np.arange(0,holdT,1), color='red', label='Uncleaned')
    plt.xlabel('Time [s]')
    plt.ylabel('Counts')
    plt.legend()
    plt.show()
#    plt.hist(data[cut]['time'] + 50, bins=np.arange(0,100,1))
#    plt.yscale('log')
#    plt.hist(data[cut]['energy']/(mass_n*grav)+(minz+1.5), bins=100)
#    plt.hist2d(data[cut]['energy']/(mass_n*grav)+(minz+1.5), bins=100)
#    plt.show()

def countHeated(data, holdT):
    cut = (data['zoff'] > 0) & (data['eStart'] < eCleanLow) & (data['energy'] > eCleanLow)
    return(np.sum(cut))

def countHeatedLost(data, holdT):
    cut = (data['zoff'] < -2.5) & (data['eStart'] < eCleanLow)
    return(np.sum(cut))

def countHeatedLostCleaner(data, holdT):
    cut = (data['zoff'] == -3) & (data['eStart'] < eCleanLow) & (data['energy'] > eCleanLow)
    return(np.sum(cut))

def countHeatedLostEscape(data, holdT):
    cut = (data['zoff'] == -4) & (data['eStart'] < eCleanLow) & (data['energy'] > eCleanLow)
    return(np.sum(cut))

def countHeatedLostHigh(data, holdT):
    cut = (data['zoff'] == -3) & (data['eStart'] < eCleanHigh) & (data['energy'] > eCleanHigh)
    return(np.sum(cut))

def countUncleanedLost(data, holdT):
    cut = (data['zoff'] == -3) & (data['eStart'] > eCleanHigh)
    return(np.sum(cut))

def countUncleaned(data, holdT):
    cleanCut = (data['zoff'] == -2) & (data['eStart'] > eCleanLow)
    bandCut = (data['eStart'] > eCleanLow)
    return(np.sum(bandCut) - np.sum(cleanCut))

def countUncleanedHigh(data, holdT):
    cleanCut = (data['zoff'] == -2) & (data['eStart'] > eCleanHigh)
    bandCut = (data['eStart'] > eCleanHigh)
    return(np.sum(bandCut) - np.sum(cleanCut))

def popBreakdown(data, holdT):
    lowCut = (data['eStart'] < eCleanLow)
    midCut = (data['eStart'] >= eCleanLow) & (data['eStart'] < eCleanHigh)
    highCut = (data['eStart'] >= eCleanHigh)
    print("LowE:", np.sum(lowCut), "MidE:", np.sum(midCut), "HighE", np.sum(highCut))

def popBreakdownDet(data, holdT):
    tStart = holdT + 20
    tEnd = holdT + 20 + 120
    lowCut = (data['eStart'] < eCleanLow) & (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEnd)
    midCut = (data['eStart'] >= eCleanLow) & (data['eStart'] < eCleanHigh) & (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEnd)
    highCut = (data['eStart'] >= eCleanHigh) & (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEnd)
    print("LowE:", np.sum(lowCut), "MidE:", np.sum(midCut), "HighE", np.sum(highCut))

def countFirstDip(data, holdT):
    tStart = holdT
    tEnd = holdT + 20
    cutUncl = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEnd) & (data['eStart'] >= eCleanLow)
    cutHeat = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEnd) & (data['eStart'] < eCleanLow)
    cutAll = (data['zoff'] > 0) & (data['time'] > tStart) & (data['time'] < tEnd)
    print("Uncl:", np.sum(cutUncl), "Heat:", np.sum(cutHeat), "All:", np.sum(cutAll))

np.random.seed(2303616184)

dt = np.dtype([('rLenFront', np.uint32, (1)),
               ('energy', np.float64, (1)),
               ('theta', np.float64, (1)),
               ('time', np.float64, (1)),
               ('settlingTime', np.float64, (1)),
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

#shortData = np.fromfile("/Volumes/SanDisk/fullTau_Short.bin", dtype=dt)
#shortData = np.fromfile("../datafiles/fulltau/fullTau_Short.bin", dtype=dt)
shortData = np.fromfile("../datafiles/fulltau/ccw-heat/fullTau_Short_ccw_clean.bin", dtype=dt)
#shortData = np.fromfile("../datafiles/fulltau/high_stat/fullTau_Short_kill.bin", dtype=dt)
#shortData = np.fromfile("../datafiles/fulltau/high_stat/cleanOnlyTau_Short_kill.bin", dtype=dt)
assert(np.all(shortData['rLenFront']==shortData[0]['rLenFront']))
#longData = np.fromfile("/Volumes/SanDisk/fullTau_Long.bin", dtype=dt)
#longData = np.fromfile("../datafiles/fulltau/fullTau_Long.bin", dtype=dt)
longData = np.fromfile("../datafiles/fulltau/ccw-heat/fullTau_Long_ccw_clean.bin", dtype=dt)
#longData = np.fromfile("../datafiles/fulltau/high_stat/fullTau_Long_kill.bin", dtype=dt)
#longData = np.fromfile("../datafiles/fulltau/high_stat/cleanOnlyTau_Long_kill.bin", dtype=dt)
assert(np.all(longData['rLenFront']==longData[0]['rLenFront']))

#shortDataNoNan = shortData[~np.isnan(shortData['energy'])]
#longDataNoNan = longData[~np.isnan(longData['energy'])]

#plotEDiff(shortData, 20)
#plotEDiff(longData, 1400)

#plotClean(shortData, 20)
#plotClean(longData, 1400)

#Uncleaned Error
#print("uncleaned   uncleanedHigh   heated   uncleanedLost   heatedLost  heatedLostCleaner  heatedLostEscape  heatedLostHigh")
#print(countUncleaned(shortData, 20), countUncleanedHigh(shortData, 20), countHeated(shortData, 20), countUncleanedLost(shortData, 20), countHeatedLost(shortData, 20), countHeatedLostCleaner(shortData, 20), countHeatedLostEscape(shortData, 20), countHeatedLostHigh(shortData, 20))
#print(countUncleaned(longData, 1400), countUncleanedHigh(longData, 1400), countHeated(longData, 1400), countUncleanedLost(longData, 1400), countHeatedLost(longData, 1400), countHeatedLostCleaner(longData, 1400), countHeatedLostEscape(longData, 1400), countHeatedLostHigh(longData, 1400))
#countFirstDip(shortData, 20)
#countFirstDip(longData, 1400)
plotDip1HighE(shortData, 20)
#plotDip1HighE(longData, 1400)
#plotDip1Heat(shortData, 20)
plotDip1Heat(longData, 1400)
#plotDip1_43cm(longData, 1400)
#plotStorageClean(shortData, 20)
#plotStorageClean(longData, 1400)


#Lifetime
numShort = countUCN(shortData, 20)
timeShort = avgT(shortData, 20)
numLong = countUCN(longData, 1400)
timeLong = avgT(longData, 1400)
#plotDip2(shortData, 20)
#plotDip2(longData, 1400)
lt = (timeLong-timeShort)/(np.log((numShort/len(shortData))/(numLong/len(longData))))
lterr = np.sqrt(1/numShort + 1/numLong)*(timeLong-timeShort)/(np.log((numShort/len(shortData))/(numLong/len(longData))))**2
print(timeShort, timeLong, numShort, numLong, (timeLong-timeShort), lt, lterr)
heatedNumber = countHeatedLostCleaner(longData, 1400)
print("num Heated:", heatedNumber)

#Systematics
#print(ratDip2(shortData, 20), ratDip2(longData, 1400))
#print(sumCleanChk(shortData, 20), sumCleanChk(longData, 1400))

#Other Checks
#popBreakdown(shortData, 20)
#popBreakdown(longData, 1400)
#popBreakdownDet(shortData, 20)
#popBreakdownDet(longData, 1400)
