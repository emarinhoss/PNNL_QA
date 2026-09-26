# MMC setup (reconnection): creates `last` folders recon_006_MR_<MI/ME>_c0_<c>.
# Each gets level 0 (ssrecon_wv_0.pin, 400x200 cells, cray_0.qsub); the first
# half also get level 1 (800x400, cray_1.qsub) and the first quarter level 2
# (1600x800, cray_2.qsub). Template: input.py. MI/ME is uniform on [25, 100].
# Next step: run_cases.sh.

import os
import random
import string

# Set SEED to an integer to make the sample reproducible
# (None seeds from the operating system, as before).
SEED = None
random.seed(SEED)
	

# get paramenters that remain the same for all runs
f = open('input.py', 'r')
info = f.read()
f.close()

# Number of runs
last = 48  # number of runs

for count in range(last):
	# ion-to-electron mass ratio MI/ME, uniform on [25, 100], with MI = 1
	v = random.uniform(25,100)
	mi = 1.0
	me = mi/v

	# uniform varying values of the speed of light
	# c = random.uniform(8,12)
	c = 1.00

	light = ("LIGHT_SPEED =" + str(c))
	masse = ("ME =" + str(me))
	massi = ("MI =" + str(mi))

	folder = ("recon_006_MR_" + str(v) + "_c0_" + str(c))


	## create Warpx .pin file
	out = open('ssrecon_wv_0.pin','w')
	## write data into .pin file
	out.write('# -*- python -*- \n')
	out.write('# The following parameters have been randomly generated. \n')
	out.write(light)
	out.write("\n")
	out.write(masse)
	out.write("\n")
	out.write(massi)
	out.write("\n")
	out.write('# -- End of randomly generated data. --')
	#out.write("rname = ssrecon_wv_0")
	out.write("\n")
	out.write("nx = 400")
	out.write("\n")
	out.write("ny = 200")
	out.write("\n")
	out.write(info)
	out.close()

	os.mkdir(folder)
	os.system("mv ssrecon_wv_0.pin " + folder)
	os.system("cp cray_0.qsub " + folder)

	if count < last//2:
		out = open('ssrecon_wv_1.pin','w')
		out.write('# -*- python -*- \n')
		out.write('# The following parameters have been randomly generated. \n')
		out.write(light)
		out.write("\n")
		out.write(masse)
		out.write("\n")
		out.write(massi)
		out.write("\n")
		out.write('# -- End of randomly generated data. --')
		#out.write("rname = ssrecon_wv_1")
		out.write("\n")
		out.write("nx = 800")
		out.write("\n")
		out.write("ny = 400")
		out.write("\n")
		out.write(info)
		out.close()
		os.system("mv ssrecon_wv_1.pin " + folder)
		os.system("cp cray_1.qsub " + folder)

	if count < last//4:
		out = open('ssrecon_wv_2.pin','w')
		out.write('# -*- python -*- \n')
		out.write('# The following parameters have been randomly generated. \n')
		out.write(light)
		out.write("\n")
		out.write(masse)
		out.write("\n")
		out.write(massi)
		out.write("\n")
		out.write('# -- End of randomly generated data. --')
		#out.write("rname = ssrecon_wv_2")
		out.write("\n")
		out.write("nx = 1600")
		out.write("\n")
		out.write("ny = 800")
		out.write("\n")
		out.write(info)
		out.close()
		os.system("mv ssrecon_wv_2.pin " + folder)
		os.system("cp cray_2.qsub " + folder)
