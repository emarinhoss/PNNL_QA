import glob
import os
import subprocess
import wxdata as wxdata2
from pylab import *
from numpy import *
import numpy



d = glob.glob("./recon_001*")

meanf = numpy.zeros(41, numpy.float)
varf  = numpy.zeros(41, numpy.float)

for m in range(0,len(d)):
	filename = os.path.join(d[m], "ssrecon_wv_0.dat")
	data = numpy.genfromtxt(filename)
	meanf = meanf + data
	varf  = varf + data*data

N = len(d)
varf = (varf/(N-1.) - meanf*meanf/(N*N-N))#/N
meanf = meanf/N

numpy.savetxt('mc_mean.dat',meanf)
numpy.savetxt('mc_vari.dat',varf)
