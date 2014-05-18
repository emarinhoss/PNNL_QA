import glob
import wxdata as wxdata2
from numpy import *

d = glob.glob("./MC_RUNS/advect_002_U_*")

meanf = zeros(1000, float)
varf  = zeros(1000, float)

Num = array([36,72,150,520,750,1150,1500,2500,3500,4500])

for k in range(len(Num)):
    mn = random.randint(4999, size=Num[k])
    
    for m in range(0,len(mn)):
        filename = os.path.join(d[mn[m]], "advect_mc")
        dh = wxdata2.WxData(filename, 10)
        q = dh.read('qnew')
        meanf += q[:,1]
        varf += q[:,1]*q[:,1]
        dh.close()
        
    N = len(d)
    varf = (varf/(N-1.) - meanf*meanf/(N*N-N))
    meanf = meanf/N
