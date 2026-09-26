# Multilevel Monte Carlo (MMC) setup: creates `last` run folders
# advect_001_U_<value>. Every folder gets a level-0 input advect_0.pin; the
# first half also get level 1 (advect_1.pin) and the first quarter level 2
# (advect_2.pin), each on a finer grid. Template: input3.py.
# The uncertain parameter is the initial ion velocity u_i, uniform on [1e-10, 1e-4].
# Also writes the sampled values to mmc_values.txt (overwriting it).
# Next step: run_cases.sh.

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
avg = open('mmc_values.txt','w')

# Number of runs
last = 36  # number of runs

for count in range(last):
	# initial ion velocity u_i, uniform on [1e-10, 1e-4]
	v = random.uniform(1.0e-10,1.0e-4)
	#v = random.gauss(3,1.0)

	ux = ("value = " + str(v))

	folder = ("advect_001_U_" + str(v))


	## create Warpx .pin file
	out = open('advect_0.pin','w')
	## write data into .pin file
	out.write('# -*- python -*- \n')
	out.write('# The following parameters have been randomly generated. \n')
	out.write(ux)
	out.write("\n")
	out.write('# -- End of randomly generated data. --')
	out.write("\n")
	#out.write("rname = advect_0")
	out.write("\n")
	out.write("nx = 100")
	out.write("\n")
	out.write("ny = 100")
	out.write("\n")
	out.write(info)
	out.close()

	os.mkdir(folder)
	os.system("mv advect_0.pin " + folder)
	avg.write(str(v))
	avg.write('\n')

	if count < last//2:
		out = open('advect_1.pin','w')
		out.write('# -*- python -*- \n')
		out.write('# The following parameters have been randomly generated. \n')
		out.write(ux)
		out.write("\n")
		out.write('# -- End of randomly generated data. --')
		out.write("\n")
		#out.write("rname = advect_1")
		out.write("\n")
		out.write("nx = 500")
		out.write("\n")
		out.write("ny = 500")
		out.write("\n")
		out.write(info)
		out.close()
		os.system("mv advect_1.pin " + folder)

	if count < last//4:
		out = open('advect_2.pin','w')
		out.write('# -*- python -*- \n')
		out.write('# The following parameters have been randomly generated. \n')
		out.write(ux)
		out.write("\n")
		out.write('# -- End of randomly generated data. --')
		out.write("\n")
		#out.write("rname = advect_2")
		out.write("\n")
		out.write("nx = 1000")
		out.write("\n")
		out.write("ny = 1000")
		out.write("\n")
		out.write(info)
		out.close()
		os.system("mv advect_2.pin " + folder)

avg.close()
