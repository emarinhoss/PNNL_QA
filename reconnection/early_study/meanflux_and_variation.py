import glob
import numpy
import os
import wxdata as wxdata2
from pylab import *
from numpy import *

frames = 10
d = glob.glob("/home/sousae/UQ_PNNL/recon_MR*")

flux  = numpy.zeros((size(d),frames+1), numpy.float)
meanf = numpy.zeros(frames+1, numpy.float)
varf  = numpy.zeros(frames+1, numpy.float)

for m in range(0,size(d)):
	filename = os.path.join(d[m], "ssrecon_wv")
	for n in range(0,frames+1):
		dh = wxdata2.WxData(filename, n)
		q = dh.read('qnew')
		dx = (q.grid.upperBounds[1]-q.grid.lowerBounds[1])/q.grid.numPhysCells[1]
		ny = int(round(q.grid.numPhysCells[1]/2))
		#ny = q.grid.numPhysCells[1]
		by = q[:, ny, 14]
		
		flux[m,n] = dx*numpy.sum( numpy.fabs(by) )
		dh.close()
		
Ly = 12.8
wce = 0.1
T = linspace(0.0, frames*2.0*wce, frames+1)
flux = flux/(2*Ly)

for l in range(0, size(d)):
	flux[l,:] = 2*wce*flux[l,:]/flux[l,0] # rescale to match GEM conditions
	meanf = meanf + flux[l,:]
	varf  = varf + flux[l,:]*flux[l,:]
	
meanf = meanf/size(d)
varf  = varf/size(d) - meanf*meanf

figure(1)
font = {'fontsize'   : 20}
#plot(T, meanf, '-b')
errorbar(T, meanf, varf, ecolor='red')
xlabel(r'$\omega_{ci}t$',font)
ylabel('Reconnected flux',font)
#legend(('TF_128x64_nocor','TF_256x128_nocor','TF_128x64_divBcor','TF_128x64_divcor','TF_512x256_divcor','HallMHD'),2)
yticks(fontsize=18)
xticks(fontsize=18)
savefig('reconflux.png')
show()
