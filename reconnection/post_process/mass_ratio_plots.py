import glob
import os
import subprocess
import wxdata as wxdata2
from pylab import *
from numpy import *
import numpy


frames = 40
Ly = 12.8
wci = 0.1
time = 400.
T = linspace(0.0, time*wci, frames+1)

########################################################################
########################################################################
# Multilevel Monte Carlo
########################################################################
########################################################################
mean_mmc = numpy.genfromtxt("./old_mmc/mean.dat")
varn_mmc = numpy.genfromtxt("./old_mmc/vari.dat")

########################################################################
########################################################################
# Monte Carlo
########################################################################
########################################################################
mean_mc = numpy.genfromtxt("./mc/mc_mean.dat")
varn_mc = numpy.genfromtxt("./mc/mc_vari.dat")

########################################################################
########################################################################
# Probabilistic Collocation Method
########################################################################
########################################################################
mean_pcm = numpy.genfromtxt("./pcm/pcm_mean_2.dat")
varn_pcm = numpy.genfromtxt("./pcm/pcm_vari_2.dat")

########################################################################
########################################################################
# Plot
########################################################################
########################################################################

figure(1)
font = {'fontsize'   : 20}
subplot(2,1,1),plot(T,mean_mmc,'c',T,mean_pcm,'r',T,mean_mc,'g',linewidth=2)
legend(('MMC','PCM','MC'),2)
#xlabel('Time ($\omega_{ci}t$)',font)
ylabel('Flux',font)
yticks(fontsize=18)
xticks(fontsize=18)


#subplot(2,1,2),plot(T,varn_mmc[:,0],'c',T,varn_mmc[:,1],'.c',T,varn_pcm,'r',T,varn_mc,'g',linewidth=2)
subplot(2,1,2),plot(T,varn_mmc[:,0],'c',T,varn_mmc[:,1],'.c',T,varn_mmc[:,2],'--c',T,varn_pcm,'r',T,varn_mc,'g',linewidth=2)
legend(('MMC_giles','MMC_mod','MMC_def','PCM','MC'),2)
xlabel('Time ($\omega_{ci}t$)',font)
ylabel('$\sigma^2$',font)
yticks(fontsize=18)
xticks(fontsize=18)

savefig('mass_ratio_recon_comparison.png')
