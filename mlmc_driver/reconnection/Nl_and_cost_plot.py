import os
import glob
import math
import numpy as np
from pylab import *
import original_mmc_warpx as mmc
import wxdata as wx

#os.system("rm -rf L*")

N = 5							# inital number of samples
M = 4
tol = [1.0e-6,0.66e-6,0.33e-6,0.1e-6]	# tolerance
nx = 25						# minimum spacial resolution
res= 1600						# plotting resolution
frames = 10

mmccost = np.zeros(len(tol), np.float)
mccost  = np.zeros(len(tol), np.float)
Nll = []

for n in range(len(tol)):
	e  = tol[n]
	ml = []
	ui = []
	suml1 = []
	suml2 = []
	suml3 = []
	suml4 = []
	suml5 = []
	L = 0

	while 1:
		# Step 1: run inital cases
		# set up runs
		ml.append(M**L*nx)
		if not(os.path.exists("L"+str(L))):
			os.mkdir("L"+str(L))
		sums = mmc.preprocess(ml,N,L,frames,res,ui)
		suml1.append(N)
		suml2.append(sums[0,:])
		suml3.append(sums[1,:])
		suml4.append(sums[2,:])
		suml5.append(sums[3,:])
	
		# Step 2: estimate the variance using the inital number of samples
		Vl = mmc.variance(suml3,suml2,suml1)
	
		# Step 3: Calculate Optimal Nl, l=0,1,...,L
		Nl = np.ceil(2*e**(-2)*np.sqrt(Vl/ml)*np.sum(np.sqrt(Vl/ml)))
	
		# Step 4: Evaluate extra samples at each level as needed for the
		# new Nl
		for k in range(0,L+1):
			Nm = Nl[k] - suml1[k]
			if Nm > 0.:
				sums = mmc.preprocess(ml,int(Nm),k,frames,res,ui)
				suml1[k] = suml1[k] + int(Nm)
				suml2[k] = suml2[k] + sums[0,:]
				suml3[k] = suml3[k] + sums[1,:]
				suml4[k] = suml4[k] + sums[2,:]
				suml5[k] = suml5[k] + sums[3,:]
	
		# Step 5: If L>=1 test for convergence using YL = M^{\alpha}
		if L>=1:
			YL1 = 1./M*suml2[L-1]/suml1[L-1]
			yl1 = max(np.abs(YL1))
			YL  = suml2[L]/suml1[L]
			yl  = max(np.abs(YL))
			converged = (max(yl1,yl)<1/(np.sqrt(2))*(M-1)*e)
			if converged:
				break
	
		# Step 6: If not converged, set L=L+1 and go back to 2
		L = L + 1
	Nll.append(Nl)
	for m in range(0,L+1):
		mmccost[n] = mmccost[n] + (Nl[m]*M**m)
		var = suml3[m]/suml1[m] - (1/(suml1[m]**2.-suml1[m]))*(suml2[m])**2
		mccost[n]  = mccost[n] + (2*var[0]/e**2)*M**m

mmccost = (1.+1./M)*mmccost


figure(1)
font = {'fontsize'   : 20}
subplot(2,1,1),semilogy(range(0,len(Nll[0])),Nll[0],'-sb',range(0,len(Nll[1])),Nll[1],'-ok',range(0,len(Nll[2])),Nll[2],'-^g',range(0,len(Nll[3])),Nll[3],'-*r',linewidth=2)
legend(('$\epsilon=1.0e-6$','$\epsilon=6.6e-7$','$\epsilon=3.3e-7$','$\epsilon=1.0e-7$'),1)
xlabel('level $\ell$',font)
ylabel('$N_l$',font)
yticks(fontsize=18)
xticks(fontsize=18)

subplot(2,1,2),loglog(tol,tol*mmccost,'-sb',tol,tol*mccost,'-ok',linewidth=2)
legend(('MMC','MC'),1)
xlabel('$\epsilon$',font)
ylabel('$\epsilon^2 cost$',font)
yticks(fontsize=18)
xticks(fontsize=18)

savefig('Nl_and_cost_plot_25_M4.png')

show()
