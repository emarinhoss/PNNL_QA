import os

# get paramenters that remain the same for all runs
f = open('input3.py', 'r')
info = f.read()
f.close

# Number samples
last = 65

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

a = 1.0e-10
b = 1.0e-4

for l in range(last):
    v = a + (b-a)*(Ps[l]+1)*0.5
    
    ux = ("value =" + str(v))

    folder = ("advect_003_U_" + str(v))


    ## create Warpx .pin file
    out = open('advect_pcm.pin','w')
    ## write data into .pin file
    out.write('# -*- python -*- \n')
    out.write('# The following parameters has been randomly generated. \n')
    out.write(ux)
    out.write("\n")
    out.write('# -- End of randomly generated data. --')
    out.write("\n")
    out.write("nx = 100")
    out.write("\n")
    out.write("ny = 100")
    out.write("\n")
    out.write(info)
    out.close

    out = open('weight.w','w')
    out.write(str(Ws[l]))
    out.close

    os.mkdir(folder)
    os.system("mv advect_pcm.pin " + folder)
    os.system("mv weight.w " + folder)
