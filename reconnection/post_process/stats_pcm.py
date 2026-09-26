# PCM statistics: reads recon_pcm.dat and weight.w in every run folder and
# writes the quadrature mean (pcm_mean_2.dat) and variance (pcm_vari_2.dat).
# Usage: python stats_pcm.py [glob]   (default glob: ./recon_002*)

import glob
import sys
import os
import subprocess
import wxdata as wxdata2
from pylab import *
from numpy import *
import numpy


frames = 40
Ly = 12.8
wci = 0.1

d = glob.glob(sys.argv[1] if len(sys.argv) > 1 else "./recon_002*")

flux  = numpy.zeros((len(d),frames+1), numpy.float)
meanf = numpy.zeros(frames+1, numpy.float)
varf  = numpy.zeros(frames+1, numpy.float)
weight = []

for m in range(0,len(d)):
        filename = os.path.join(d[m], "recon_pcm.dat")
        print d[m]
        weigthfile=os.path.join(d[m], "weight.w")
        wts = open(weigthfile,'r')
        weight.append(float(wts.readline()))
        wts.close()
        data = genfromtxt(filename)
        meanf = meanf + data*weight[m]
        varf = varf + data*data*weight[m]

varf = varf - meanf*meanf
        
numpy.savetxt('pcm_mean_2.dat',meanf)
numpy.savetxt('pcm_vari_2.dat',varf)
