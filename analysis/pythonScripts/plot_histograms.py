import sys

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

jtonev = 6.2415091e27

if(len(sys.argv) != 2):
	sys.exit("Error! Usage: ./plotArrivalTime fname")

fname = sys.argv[1]

dt = np.dtype([('buf', np.int32, (1)),
			   ('energy', np.float64, (1)),
               ('theta', np.float64, (1)),
               ('time', np.float32, (50)),
               ('nhit', np.int32, (50)),
               ('buf2',np.int32,(1))])

data = np.fromfile(fname,dtype=dt) 

energy= []
theta = []
meanT = []
stdT  = []
energyHist = []
thetaHist = []
tDiffHist  = []

# find maximum and minimum energy for hist bins
minE = np.floor(min(data['energy'])*jtonev)
maxE = np.ceil(max(data['energy'])*jtonev)
#print(minE,maxE)

for t in data:
	
	theta.append(t['theta'])
	eTmp = t['energy']*jtonev
	tInt = []
	tPrev = 0.0
	for i,tTmp in enumerate(t['time']):
		if t['nhit'][i] == 2:# or t['nhit'][i] == 1:
			if not tPrev == 0.0:
				energyHist.append(eTmp)
				thetaHist.append(t['theta'])
				tDiffHist.append(tTmp - tPrev)
				tInt.append(tTmp - tPrev)
			tPrev = tTmp
			
	#if eTmp > 38:
	#	continue
	energy.append(t['energy'] * (jtonev))
	if len(tInt) > 0:
		#meanT.append(np.mean(tInt))
		meanT.append(max(tInt))
		stdT.append(np.std(tInt))
	else:
		meanT.append(np.max(t['time']))
		#meanT.append(-1.0)
		stdT.append(0.0)

xedgesEn = np.linspace(float(minE),float(maxE),(maxE-minE))
if max(tDiffHist) < 2.0:
	yedgesT = np.linspace(0.0,2.0,50)
else:
	yedgesT = np.linspace(0.0,float(np.ceil(max(tDiffHist))),50)
Hist, xedgesEn,yedgesT = np.histogram2d(energyHist, tDiffHist, bins=(xedgesEn,yedgesT))
Hist = Hist.T # Fixes a weird mismatch error
XEn, YT = np.meshgrid(xedgesEn,yedgesT)
plt.figure(1)
plt.pcolormesh(XEn,YT,Hist,norm=colors.LogNorm(vmin=1,vmax=Hist.max()))
plt.colorbar()
plt.title('Time to Next TD Hit Histogram')
plt.xlabel('Energy (neV)')
plt.ylabel('Time to next bounce (s)')

plt.figure(2)
plt.errorbar(energy,meanT,yerr=stdT, fmt = 'b.')
plt.title('Time to Next TD Hit Raw')
plt.xlabel('Energy (neV)')
plt.ylabel('Avg. time to next bounce (s)')

xedgesA = np.linspace(0.0,np.pi / 2,50)
if max(tDiffHist) < 2.0:
	yedgesT = np.linspace(0.0,2.0,50)
else:
	yedgesT = np.linspace(0.0,float(np.ceil(max(tDiffHist))),50)
Hist, xedgesA,yedgesT = np.histogram2d(thetaHist, tDiffHist, bins=(xedgesA,yedgesT))
Hist = Hist.T # Fixes a weird mismatch error
XA, YT = np.meshgrid(xedgesA,yedgesT)
plt.figure(3)
plt.pcolormesh(XA,YT,Hist,norm=colors.LogNorm(vmin=1,vmax=Hist.max()))
plt.colorbar()
plt.title('Time to Next TD Hit Angular Histogram')
plt.xlabel('Initial angle ')
plt.ylabel('Time to next bounce (s)')


fig = plt.figure(4)
ax = Axes3D(fig)#fig.add_subplot(111,projection='3d')
#plt.axes(projection="3d")

#ax.scatter(energyHist,thetaHist,tDiffHist,c=tDiffHist)
ax.scatter(energy, theta, meanT, c=meanT)
plt.xlabel('Initial Energy (neV)')
plt.ylabel('Initial Angle')
plt.title('Time to Next TD Hit Scatter Plot')
plt.show()


#plt.show()


#print(np.mean(meanT), np.std(meanT))
#plt.figure(1)
#plt.errorbar(energy,meanT,yerr=stdT, fmt = 'b.')


#plt.ylabel('Avg. time to next bounce (s)')


