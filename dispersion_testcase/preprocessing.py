# Monte Carlo (MC) setup: creates `last` run folders advect_002_U_<value>,
# each holding advect_mc.pin = sampled value + grid size + the input3.py template.
# The uncertain parameter is the initial ion velocity u_i, uniform on [1e-10, 1e-4].
# Also writes the sampled values to mc_values.txt (overwriting it).
# Next step: run_mc.sh.

import os
import random
import string

# Set SEED to an integer to make the sample reproducible
# (None seeds from the operating system, as before).
SEED = None
random.seed(SEED)
	

# get paramenters that remain the same for all runs
f = open('input3.py', 'r')
info = f.read()
f.close()

#write all the values on a file
avg = open('mc_values.txt','w')

# Number of runs
last = 63  # number of runs

for count in range(last):
	# initial ion velocity u_i, uniform on [1e-10, 1e-4]
	v = random.uniform(1.0e-10,1.0e-4)
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
	avg.write(str(v))
	avg.write('\n')

avg.close()
