# MC statistics: reads ssrecon_wv_0.dat (from calc_flux.py) in every run
# folder and writes the sample mean (mc_mean.dat) and unbiased sample
# variance (mc_vari.dat) of the reconnected flux per frame.
# Usage: python stats_mc.py [glob]   (default glob: ./recon_001*)

import glob
import sys
import os
import subprocess
import wxdata as wxdata2
from pylab import *
from numpy import *
import numpy



d = glob.glob(sys.argv[1] if len(sys.argv) > 1 else "./recon_001*")

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
