def preprocess(m,N,L,frames,res):
	import numpy as np
	import os
	import random
	import fileinput
	import sys
	import string
	import wxdata as wx
	
	x = np.linspace(0.0, 1.0, res)
	sums = np.zeros((4,res),np.float)

	for count in range(int(N)):
		# uniform varying values
		v = random.uniform(1.0e-10,1.0e-4)
		#v = random.gauss(3,1.0)
		os.system("cp input3.inp mmc_testcase.inp")
		for line in fileinput.input("mmc_testcase.inp", inplace=1):
			if "        u_i = value" in line:
				line = line.replace("        u_i = value","        u_i = "+str(v))
			elif "      Cells = [nx]" in line:
				line = line.replace("      Cells = [nx]","      Cells = ["+str(m[L])+"]")
			sys.stdout.write(line)
		
		os.system("mpirun -n 1 $warpxp -i mmc_testcase.inp")
		
		df = wx.WxData('mmc_testcase', frames)
		qf = df.read('qnew')
		vf = qf[:, 1]/qf[:, 0]
		xf = np.linspace(qf.grid.lowerBounds[0],qf.grid.upperBounds[0],qf.grid.numPhysCells[0])
		df.close
		Xf = np.interp(x,xf,vf)
		
		folder = ("U_" + str(v))
		where1 = ('L'+str(L) + '/' + folder)
		os.mkdir(where1)
		os.system("mv *.h5 *.log mmc_testcase.inp " + where1)
		
		if L > 0:
			os.system("cp input3.inp mmc_testcase.inp")
			for line in fileinput.input("mmc_testcase.inp", inplace=1):
				if "        u_i = value" in line:
					line = line.replace("        u_i = value","        u_i = "+str(v))
				elif "      Cells = [nx]" in line:
					line = line.replace("      Cells = [nx]","      Cells = ["+str(m[L-1])+"]")
				sys.stdout.write(line)
				
			os.system("mpirun -n 1 $warpxp -i mmc_testcase.inp")
			
			dc = wx.WxData('mmc_testcase', frames)
			qc = dc.read('qnew')
			vc = qc[:, 1]/qc[:, 0]
			xc = np.linspace(qc.grid.lowerBounds[0],qc.grid.upperBounds[0],qc.grid.numPhysCells[0])
			dc.close
			Xc = np.interp(x,xc,vc)
			where2 = ('L'+str(L-1)+'/'+folder)
			os.mkdir(where2)
			os.system("mv *.h5 *.log mmc_testcase.inp " + where2)
		else:
			Xc=np.zeros(res, np.float)

		sums[0,:] = sums[0,:] + (Xf-Xc)
		sums[1,:] = sums[1,:] + (Xf-Xc)*(Xf-Xc)
		sums[2,:] = sums[2,:] + Xf
		sums[3,:] = sums[3,:] + Xf*Xf
	
	return sums
	
def variance(sum1,sum2,N):
	import numpy as np

	var = np.zeros(len(N),np.float)
	for k in range(0,len(N)):
		v = sum1[k]/N[k]
		m = sum2[k]/N[k]
		var[k] = max(v - m*m)
	return var
