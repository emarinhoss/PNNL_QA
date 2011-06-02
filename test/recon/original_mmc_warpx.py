import numpy as np
import os
import random
import fileinput
import sys
import string
import wxdata as wx
import subprocess

def preprocess(m,N,L,frames,ui):

	sums = np.zeros((4,frames+1),np.float)
	#home = os.getcwd()

	for count in range(1,N+1):
		# uniform varying values
		v = random.uniform(25,100)
		#v = random.gauss(3,1.0)
		ui.append(v)
		
		Xf = getdata(v,1.0,0.1,m[L],L,frames)
		
		if L > 0:
			Xc = getdata(v,1.0,0.1,m[L-1],L-1,frames)
		else:
			Xc=np.zeros(frames+1, np.float)

		sums[0,:] = sums[0,:] + (Xf-Xc)
		sums[1,:] = sums[1,:] + (Xf-Xc)*(Xf-Xc)
		sums[2,:] = sums[2,:] + Xf
		sums[3,:] = sums[3,:] + Xf*Xf
	
	return sums
	
def variance(sum1,sum2,N):
	import numpy as np

	var = np.zeros(len(N),np.float)
	for k in range(0,len(N)):
		v = sum1[k]/(N[k]-1)
		m = (sum2[k])**2/(N[k]**2-N[k])
		var[k] = max(v - m*m)
	return var

def getdata(v,LS,BFP,res,L,frames):

	os.system("cp input_mmc.inp mmc_testcase.inp")
	resx = res
	resy = res/2
	for line in fileinput.input("mmc_testcase.inp", inplace=1):
		if "ME" in line:
			line = line.replace("ME",str(v))
		elif "      Cells = [resx, resy]" in line:
			line = line.replace("      Cells = [resx, resy]","      Cells = ["+str(resx)+", "+str(resy)+"]")
		elif "LS" in line:
			line = line.replace("LS",str(LS))
		elif "bpert" in line:
			line = line.replace("bpert",str(BFP))
		sys.stdout.write(line)
		
	folder = ("U_" + str(v))
	where1 = ('L'+str(L) + '/' + folder)
	os.mkdir(where1)

	p = subprocess.Popen("crun -nodes 1 -cores 1 $warpx -i mmc_testcase.inp", shell=True)
	sts = os.waitpid(p.pid, 0)[1]
	
	fluxtfh = np.zeros(frames+1, np.float)
	
	for i in range(0,frames+1):
		# open file for reading
		dtfh = wx.WxData('mmc_testcase', i)
		qtfh = dtfh.read('qnew')        
		# get Y-axis spacing
		dxtfh = qtfh.grid.dx[0]
		nytfh = qtfh.grid.numPhysCells[1]
		bytfh = qtfh[:, int(np.ceil(nytfh/2)), 14]
		fluxtfh[i] = dxtfh*np.sum( np.fabs(bytfh) )
		dtfh.close()
	
	Ly = 12.8	
	fluxtfh = fluxtfh/(2*Ly)
	fluxtfh = 0.2*fluxtfh/fluxtfh[0] # rescale to match GEM conditions
	
	os.system('mv *.h5 *.log mmc_testcase.inp ' + where1)
	return fluxtfh
