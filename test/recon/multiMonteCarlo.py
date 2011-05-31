import os
import glob
import math
import numpy as np
from pylab import *
import original_mmc_warpx as mmc
import wxdata as wx
import time

os.system("rm -rf L*")

N = 5		# inital number of samples
M = 2		# 
e = 1.0e-2	# tolerance
L = 0		# inital level
nx = 36		# spacial resolution
frames = 10

ml = []
ui = []
suml1 = []
suml2 = []
suml3 = []
suml4 = []
suml5 = []

while 1:
	# Step 1: run inital cases
	# set up runs
	ml.append(M**L*nx)
	if not(os.path.exists("L"+str(L))):
		os.mkdir("L"+str(L))
	sums = mmc.preprocess(ml,N,L,frames,ui)
	suml1.append(N)
	suml2.append(sums[0,:])
	suml3.append(sums[1,:])
	suml4.append(sums[2,:])
	suml5.append(sums[3,:])
	
	# Step 2: estimate the variance using the inital number of samples
	Vl = mmc.variance(suml3,suml2,suml1)
	
	# Step 3: Calculate Optimal Nl, l=0,1,...,L
	Nl = np.ceil(2*e**(-2)*np.sqrt(Vl/ml)*np.sum(np.sqrt(Vl/ml)))
	print Nl
	time.sleep(5)
	
	# Step 4: Evaluate extra samples at each level as needed for the
	# new Nl
	for k in range(0,L+1):
		Nm = Nl[k] - suml1[k]
		if Nm > 0.:
			sums = mmc.preprocess(ml,int(Nm),k,frames,ui)
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
		if converged or (L>5):
			break
	
	# Step 6: If not converged, set L=L+1 and go back to 2
	L = L + 1

mmc_sol = np.zeros(frames+1, np.float)
mmc_varn= np.zeros(frames+1, np.float)
T = np.linspace(0.,40.,frames+1)
levels = []

# Multi-level Monte Carlo Solution
for m in range(0,L+1):	
	mmc_sol = mmc_sol + suml2[m]/suml1[m]
	mmc_varn= mmc_varn+ 1./suml1[m]*(suml3[m]/(suml1[m]-1) - (1/(suml1[m]**2-suml1[m]))*(suml2[m])**2)
	levels.append(m)

figure(1)
font = {'fontsize'   : 20}
errorbar(T,mmc_sol,mmc_varn)
xlabel('Position',font)
ylabel('Error',font)
yticks(fontsize=18)
xticks(fontsize=18)
savefig('Solution.png')
