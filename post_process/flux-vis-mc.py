import glob
import os
from pylab import *
from numpy import *



d = glob.glob("./recon_001*")
d.sort()
time = linspace(0,40,41)
mean = zeros(41, float)
varn = zeros(41, float)

for m in range(0,len(d)):
	filename = os.path.join(d[m], "ssrecon_wv_0.dat")
        print d[m]
	data = genfromtxt(filename)
        mean = mean + data
	varn = varn + data*data
        red = float(m)/float(len(d))
        green = 0.5
        blue = float(len(d)-m)/float(len(d))
        figure(1),plot(time,data,color=(red,green,blue))

N = len(d)
varn = varn/(N-1) - mean*mean/(N*N-N)
mean = mean/N

figure(1),errorbar(time,mean,sqrt(varn),color='g',lw=2)
font = {'fontsize'   : 20}
xlabel('Time ($\omega_{ci}t$)',font)
ylabel('Flux',font)
yticks(fontsize=18)
xticks(fontsize=18)

savefig('Flux_mass_mc.png')
