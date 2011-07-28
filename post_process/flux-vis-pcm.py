import glob
import os
from pylab import *
from numpy import *



d = glob.glob("./recon_002*")
d.sort()
time = linspace(0,40,41)
mean = zeros(41, float)
varn = zeros(41, float)
weight = []

for m in range(0,len(d)):
	filename = os.path.join(d[m], "recon_pcm.dat")
        print d[m]
	weigthfile=os.path.join(d[m], "weight.w")
	wts = open(weigthfile,'r')
	weight.append(float(wts.readline()))
	wts.close()
	data = genfromtxt(filename)
        mean = mean + data*weight[m]
	varn = varn + data*data*weight[m]
        red = float(m)/float(len(d))
        green = 0.5
        blue = float(len(d)-m)/float(len(d))
        figure(1),plot(time,data,color=(red,green,blue))


figure(1),errorbar(time,mean,sqrt(varn),mfc='red',mec='green')
font = {'fontsize'   : 20}
xlabel('Time ($\omega_{ci}t$)',font)
ylabel('Flux',font)
yticks(fontsize=18)
xticks(fontsize=18)

savefig('Flux_mass_pcm_002.png')
