# Monte Carlo (MC) setup: creates `last` run folders advect_002_U_<value>,
# each holding advect_mc.pin = sampled value + grid size + the input2.py template.
# The uncertain parameter is the amplitude a of the initial condition a*sin(x), uniform on [0.5, 3].
# Next step: run_mc.sh.

import os
import random
import string

# Set SEED to an integer to make the sample reproducible
# (None seeds from the operating system, as before).
SEED = None
random.seed(SEED)
	

# get paramenters that remain the same for all runs
f = open('input2.py', 'r')
info = f.read()
f.close()

# Number of runs
last = 63  # number of runs

for count in range(last):
	# amplitude of the initial sine wave, uniform on [0.5, 3]
	v = random.uniform(0.5,3)
	#v = random.gauss(3,1.0)

	ux = ("value =" + str(v))

	folder = ("advect_002_U_" + str(v))


	## create Warpx .pin file
	out = open('advect_mc.pin','w')
	## write data into .pin file
	out.write('# -*- python -*- \n')
	out.write('# The following parameters have been randomly generated. \n')
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
	out.close()

	os.mkdir(folder)
	os.system("mv advect_mc.pin " + folder)
