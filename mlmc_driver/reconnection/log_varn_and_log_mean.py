import os
import math
import numpy as np
from pylab import *
import original_mmc_warpx as mmc

M = 4
res = 25
L = 4
frames = 10

Pl_mean = np.zeros((L+1,res), np.float)
Pl1_mean= np.zeros((L+1,res), np.float)
Pl_sigm = np.zeros((L+1,res), np.float)
Pl1_sigm= np.zeros((L+1,res), np.float)

levels = []
uii = []
mll = []
NN = 50

for m in range(0,L+1):
	mll.append(M**m*res)
	sums = mmc.preprocess(mll,NN,m,frames,res,uii)
	Pl_mean[m,:]  = sums[2,:]/NN
	Pl_sigm[m,:]  = sums[3,:]/(NN-1) - (1/(NN**2-NN))*(sums[2,:])**2
	Pl1_mean[m,:] = sums[0,:]/NN
	Pl1_sigm[m,:] = sums[1,:]/(NN-1) - (1/(NN**2-NN))*(sums[0,:])**2
	levels.append(m)

figure(2)
font = {'fontsize'   : 20}
subplot(2,1,1),plot(levels,np.log(np.abs(Pl_mean[:,0]))/np.log(M),'-s',levels[1:L+1],np.log(np.abs(Pl1_mean[1:L+1,0]))/np.log(M),'-*r',linewidth=2)
legend(('$P_l$','$P_l-P_{l-1}$'),3)
xlabel('$\ell$',font)
ylabel('$log_M |mean|$',font)

font = {'fontsize'   : 20}
subplot(2,1,2),plot(levels,np.log(Pl_sigm[:,0])/np.log(M),'-s',levels[1:L+1],np.log(Pl1_sigm[1:L+1,0])/np.log(M),'-*r',linewidth=2)
legend(('$P_l$','$P_l-P_{l-1}$'),3)
xlabel('$\ell$',font)
ylabel('$log_M variance$',font)

savefig('log_mean_and_log_variance_plots_25_M4.png')
show()
