import os
import random
import string
	

# get paramenters that remain the same for all runs
f = open('input2.py', 'r')
info = f.read()
f.close

# Number of runs
last = 63  # number of runs

for count in range(last):
	# uniform varying values of mass ration from
	# a to b given a certain value of ME
	v = random.uniform(0.5,3)
	#v = random.gauss(3,1.0)

	ux = ("value =" + str(v))

	folder = ("advect_002_U_" + str(v))


	## create Warpx .pin file
	out = open('advect_mc.pin','w')
	## write data into .pin file
	out.write('# -*- python -*- \n')
	out.write('# The following parameters has been randomly generated. \n')
	out.write(ux)
	out.write("\n")
	out.write('# -- End of randomly generated data. --')
	out.write("\n")
	#out.write("rname = advect_mc")
	out.write("\n")
	out.write("nx = 100")
	out.write("\n")
	out.write("ny = 100")
	out.write("\n")
	out.write(info)
	out.close

	os.mkdir(folder)
	os.system("mv advect_mc.pin " + folder)
