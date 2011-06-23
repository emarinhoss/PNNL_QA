import os

# get paramenters that remain the same for all runs
f = open('input2.py', 'r')
info = f.read()
f.close

# Number samples
last = 33

pts = open(("points"+str(last)),'r')
wts = open(("weights"+str(last)),'r')

Ps=[]
Ws=[]

for k in range(last):
    p = pts.readline()
    Ps.append(float(p))

    w = wts.readline()
    Ws.append(float(w))

pts.close
wts.close

MI= 1.0
a = 25.0
b = 100.0

for l in range(last):
    v = a + (b-a)*(Ps[l]+1)*0.5
    v = MI/v
    ux = ("value =" + str(v))

    folder = ("recon_001_ME_" + str(v))


    ## create Warpx .pin file
    out = open('recon_pcm.pin','w')
    ## write data into .pin file
    out.write('# -*- python -*- \n')
    out.write('# The following parameters has been randomly generated. \n')
    out.write(ux)
    out.write("\n")
    out.write('# -- End of randomly generated data. --')
    out.write("\n")
    out.write(info)
    out.close

    out = open('weight.w','w')
    out.write(str(Ws[l]))
    out.close

    os.mkdir(folder)
    os.system("mv recon_pcm.pin weight.w " + folder)
    os.system("cp batch_pcm.msub " + folder)
