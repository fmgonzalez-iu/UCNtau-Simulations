#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.optimize import curve_fit
import csv
import sys

# Fitting function single
def sexp(x, a, b):
	return a * np.exp(-x/b)
	
guessS = (1.0,880.0)

# Fitting function double
def dexp(x, a, b, c):
	return a * np.exp(-x/b) + c * np.exp(-x/877.7)
	
guessD = (0.5,877.7,0.5)

def ltShift(data):
    fracLost = num_lost/population_size #Fraction of UCN lost during storage (with no beta decay)
    fracLostError = np.sqrt(num_lost)/population_size #Assume statistics dominated by # lost
    shift = 1380.0/np.log(1/(np.exp(-1380.0/877.7)*(1-fracLost)))-877.7
    shiftErr = fracLostError*(1380.0/(1-fracLost))/(np.log(1/(np.exp(-1380.0/877.7)*(1-fracLost)))**2)
    return (shift, shiftErr)


if(len(sys.argv) < 3):
	sys.exit("Error! Usage: python plot_norm_yields.py fNameMC fNameData")
	
time = []
nyield = []
nerror = []
g1yield = []
g2yield = []
g3yield = []
g1error = []
g2error = []
g3error = []

with open(sys.argv[1]) as csvMCFile:
	yreader = csv.reader(csvMCFile)
	for row in yreader:
		time.append(float(row[0]))
		nyield.append(float(row[1]))
		nerror.append(float(row[2]))
		g1yield.append(float(row[3]))
		g1error.append(float(row[4]))
		g2yield.append(float(row[5]))
		g2error.append(float(row[6]))
		g3yield.append(float(row[7]))
		g3error.append(float(row[8]))

timeDat = []
nyDat = []
neDat = []

with open(sys.argv[2]) as csvDataFile:
	yreader = csv.reader(csvDataFile)
	for row in yreader:
		timeDat.append(float(row[0]))
		nyDat.append(float(row[1]))
		neDat.append(float(row[2]))

# Hardcoded in times.
timeAvg = [20.0, 50.0, 100.0, 200.0, 600.0, 1200.0, 1550.0, 2000.0, 2500.0]
nyAvg = []
neAvg = []

# Normalize data -- start by summing times
for j in range(0,len(timeAvg)):
	nCounts = 0.0
	nyAvg.append(0.0)
	neAvg.append(0.0)
	for i in range(0,len(nyDat)):
		if (np.abs(timeDat[i] - timeAvg[j]) / timeDat[i] < 0.1):
			#print ((timeDat[i] - timeAvg[j]) / timeAvg[j])
			nyAvg[j] += nyDat[i]
			neAvg[j] += neDat[i]*neDat[i]
			#print timeDat[i],nyAvg[j],neAvg[j]
			nCounts = nCounts + 1.0
	neAvg[j] = np.sqrt(neAvg[j])/nCounts
	nyAvg[j] = nyAvg[j]/nCounts
	
	#print timeAvg[j],nyAvg[j], neAvg[j]
		
# Scale to same value as before
norm = nyAvg[0]
for i in range(0,len(timeAvg)):		
	nyAvg[i] = nyAvg[i] / norm
	neAvg[i] = neAvg[i] / norm

norm = nyield[0]
for i in range(0,len(time)):		
	nyield[i] = nyield[i] / norm# * np.exp(-time[i]/877.7)/np.exp(-time[0]/877.7)
	nerror[i] = nerror[i] / norm# * np.exp(-time[i]/877.7)/np.exp(-time[0]/877.7)
	


#time = (20.0, 50.0, 100.0, 200.0, 600.0, 1550.0)
##Monte Carlo Data
##E^1.30448, Emin=0.076071, cos^1.264079, thetamax=1.3, p(upscatter)=1.0, B-taper=0.0160763

#nyield = (1.0, 0.9431665, 0.8576234, 0.7148234, 0.3755098, 0.1040005)
#nerror = (0.002017432, 0.001952717, 0.001853756, 0.001674789, 0.001176051, 0.0005940558)

#g1yield = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
#g1error = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

#g2yield = (0.03609346, 0.03129722, 0.02689870, 0.01960239, 0.006105438, 0.0006022147)
#g2error = (4.433451e-4, 4.119845e-4, 3.826339e-4, 3.225143e-4, 1.798000e-4, 0.5458797e-4)

#g3yield = (0.9639065, 0.9118693, 0.8307247, 0.6952211, 0.3694043, 0.1033983)
#g3error = (1.968115e-3, 1.908762e-3, 1.813936e-3, 1.642851e-3, 1.162226e-3, 0.5915424e-3)

#E^1.30448, Emin=0.076071, cos^1.264079, thetamax=0.5, p(upscatter)=1.0, B-taper=0.0160763
#nyield = (1.0, 0.9253309, 0.8239347, 0.6589763, 0.3060428, 0.07348011)
#nerror = (0.003600714, 0.003442710, 0.003222298, 0.002834169, 0.001801626, 0.0007923003)

#g1yield = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
#g1error = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

#g2yield = (0.07088601, 0.06234399, 0.05188001, 0.03798302, 0.01178169, 0.001122338)
#g2error = (1.125212e-3, 1.058632e-3, 9.603271e-4, 8.181111e-4, 4.462559e-4, 1.281876e-4)

#g3yield = (0.9291140, 0.8629870, 0.7720547, 0.6209932, 0.2942611, 0.07235778)
#g3error = (3.420386e-3, 3.275904e-3, 3.075870e-3, 2.713523e-3, 1.745483e-3, 0.7818617e-3)

# Fit curve (to single/double)
optS,covS = curve_fit(sexp,time,nyield,p0=guessS,sigma=nerror,absolute_sigma=True)
#optS,covS = curve_fit(sexp,timeAvg,nyAvg,p0=guessS,sigma=neAvg,absolute_sigma=False)
errS = np.sqrt(np.diag(covS))

#optD,covD = curve_fit(dexp,time,nyield,p0=guessD,sigma=nerror)
#errD = np.sqrt(np.diag(covD))

for i in range(len(optS)):
	print(str(optS[i])+ ' +- ' +str(errS[i]))

#for i in range(len(optD)):
#	print(str(optD[i])+ ' +- ' +str(errD[i]))
	
timefull = np.linspace(10.0,1600.0,1000)
fitfcnS = sexp(timefull,*optS)
#fitfcnD = dexp(timefull,*optD)

plt.figure()

plt.semilogy(timefull,fitfcnS,'r',label=('$e^{-t/ \\tau }$, $ \\delta \\tau = $ ' + ('%.2f' % (optS[1]-877.7)) + ' $\pm$ ' + ('%.2f' % errS[1]) + ' s'))
#lt.plot(timefull,fitfcnD,'g',label=('$A e^{-t/ \\tau } + B e^{-t/877.7}$, $ \\tau = $ ' + ('%.2f' % (optD[1])) + ' $\pm$ ' + ('%.2f' % (errD[1])) + ' s'))
plt.errorbar(time,nyield,yerr=nerror,fmt = '.',color ='b', label='MC Yields')
#plt.errorbar(timeAvg,nyAvg,yerr=neAvg,fmt = '.',color ='orange', label='Data')

plt.title('Neutron yield with fully reflecting block')
plt.legend(loc = 'lower left')
plt.xlabel('Holding time (s)')
plt.ylabel('Neutron yield (arb. units, normalized to 20s hold)')
plt.axis([0,1600,5e-2,1.1])

plt.show()
