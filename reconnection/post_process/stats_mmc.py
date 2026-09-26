import glob
import os
import subprocess
import wxdata as wxdata2
from numpy import *
from pylab import *
import numpy

dir  = glob.glob("./recon_001*")
time = linspace(0.0,40.0,41)

mean = numpy.zeros((3,41), numpy.float)
varn = numpy.zeros((3,41), numpy.float)
var2 = numpy.zeros((3,41), numpy.float)
N = numpy.zeros(3, numpy.float)

aver = numpy.zeros(41, numpy.float)
ave2 = numpy.zeros(41, numpy.float)

vari = numpy.zeros(41, numpy.float)
varg = numpy.zeros(41, numpy.float)

for m in range(0,len(dir)):

    if os.path.exists(os.path.join(dir[m],'ssrecon_wv_2.pin')):
        level2 = numpy.genfromtxt(os.path.join(dir[m], "ssrecon_wv_2.dat"))
        level1 = numpy.genfromtxt(os.path.join(dir[m], "ssrecon_wv_1.dat"))
        level0 = numpy.genfromtxt(os.path.join(dir[m], "ssrecon_wv_0.dat"))
        flux2 = level2-level1
        flux1 = level1-level0
        flux0 = level0
        mean[2,:] = mean[2,:] + flux2
        mean[1,:] = mean[1,:] + flux1
        mean[0,:] = mean[0,:] + flux0
        varn[2,:] = varn[2,:] + flux2*flux2
        varn[1,:] = varn[1,:] + flux1*flux1
        varn[0,:] = varn[0,:] + flux0*flux0
        var2[2,:] = var2[2,:] + level2*level2 - level1*level1
        var2[1,:] = var2[1,:] + level1*level1 - level0*level0
        var2[0,:] = var2[0,:] + level0*level0
#        var2[2,:] = var2[2,:] + flux2*flux2 - flux1*flux1
#        var2[1,:] = var2[1,:] + flux1*flux1 - flux0*flux0
#        var2[0,:] = var2[0,:] + flux0*flux0
        N[2] = N[2] + 1.0
        N[1] = N[1] + 1.0
        N[0] = N[0] + 1.0
        plot(time,level2,'--c')
        plot(time,level1,'--k')
        plot(time,level0,'--m')
    elif os.path.exists(os.path.join(dir[m],'ssrecon_wv_1.pin')):
        level1 = numpy.genfromtxt(os.path.join(dir[m], "ssrecon_wv_1.dat"))
        level0 = numpy.genfromtxt(os.path.join(dir[m], "ssrecon_wv_0.dat"))
        flux1 = level1-level0
        flux0 = level0
        mean[1,:] = mean[1,:] + flux1
        mean[0,:] = mean[0,:] + flux0
        varn[1,:] = varn[1,:] + flux1*flux1
        varn[0,:] = varn[0,:] + flux0*flux0
        var2[1,:] = var2[1,:] + level1*level1 - level0*level0
        var2[0,:] = var2[0,:] + level0*level0
#        var2[1,:] = var2[1,:] + flux1*flux1 - flux0*flux0
#        var2[0,:] = var2[0,:] + flux0*flux0
        N[1] = N[1] + 1.0
        N[0] = N[0] + 1.0
        plot(time,level1,'--k')
        plot(time,level0,'--m')
    else:
	level0 = numpy.genfromtxt(os.path.join(dir[m], "ssrecon_wv_0.dat"))
        flux0 = level0
        mean[0,:] = mean[0,:] + flux0
        varn[0,:] = varn[0,:] + flux0*flux0
        var2[0,:] = var2[0,:] + level0*level0
#        var2[0,:] = var2[0,:] + flux0*flux0
        N[0] = N[0] + 1.0
        plot(time,level0,'--m')

for l in range(0,3):
    aver = aver + mean[l,:]/N[l]
    ave2 = ave2 + var2[l,:]/N[l]
    varg = varg + (varn[l,:]/(N[l]-1.0) - mean[l,:]*mean[l,:]/(N[l]*N[l]-N[l]))/N[l]
    
ave2 = ave2 - aver*aver

vari = varn[0,:]/(N[0]-1.0) - mean[0,:]*mean[0,:]/(N[0]*N[0]-N[0])
for m in range(1,3):    
    vari = vari + (varn[l,:]/(N[l]-1.0) - mean[l,:]*mean[l,:]/(N[l]*N[l]-N[l]))/N[l]
    #vari = vari + (varn[l,:] - mean[l,:]*mean[l,:]/N[l])/N[l]

errorbar(time,aver,sqrt(abs(vari)))
font = {'fontsize'   : 20}
xlabel('Time ($\omega_{ci}t$)',font)
ylabel('Flux',font)
yticks(fontsize=18)
xticks(fontsize=18)
savefig('Flux_mass_mmc.png')

numpy.savetxt('mean.dat',aver)
DataOut = column_stack((varg,vari,ave2))
numpy.savetxt('vari.dat',DataOut)
