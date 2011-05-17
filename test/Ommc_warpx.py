import numpy as np
import os
import random
import string
import wxdata as wx

def preprocess(m,N):
	n = len(m)
	L = n-1
	frames = 10
	x = np.linspace(0.0, 1.0, m[0])
	sums = np.zeros((4,m[0]),np.float)

	# get paramenters that remain the same for all runs
	f = open('mmc_testcase.inp', 'r')
	info = f.read()
	f.close

	for count in range(N):
		# uniform varying values
		v = random.uniform(0.5,3)
		#v = random.gauss(3,1.0)

		ux = ("a = " + str(v))
		
		folder = ("U_" + str(v))

		## create Warpx .pin file
		out = open('testcase.inp','w')
		for line in info:
			out.write(line.replace("a = value",ux))
		out.close

		where1 = ('L'+str(L) + '/' + folder)
		os.mkdir(where1)
		os.system("mpirun -n 1 $warpxp -i testcase.inp")
		os.system("mv *.h5 *.log " + where1)
		
		
		df = wx.WxData(where1+'/euler_mmc', frames)
		qf = dh.read('qnew')
		vf = qf[:, 1]/qf[:, 0]
		xf = np.linspace(qf.grid.lowerBounds[0],qf.grid.upperBounds[0],qf.grid.numPhysCells[0])
		df.close
		Xf = np.interp(x,xh,vl)

		
		Xc=np.zeros(m[0], np.float)
		
		sums[1,:] = sums[1,:] + np.sum(Xf-Xc)
		sums[2,:] = sums[2,:] + np.sum((Xf-Xc)*(Xf-Xc))
		sums[3,:] = sums[3,:] + np.sum(Xf)
		sums[4,:] = sums[4,:] + np.sum(Xf*Xf)
	
	return sums
