import os
import random
import string
	

# get paramenters that remain the same for all runs
f = open('input.py', 'r')
info = f.read()
f.close

# Number of runs
last = 10  # number of runs

for count in range(last):
	# uniform varying values of mass ration from
	# a to b given a certain value of ME
	# v = random.uniform(25,100)
	v = 25.0
	mi = 1.0
	me = mi/v

	# uniform varying values of the speed of light
	c = random.uniform(1.0,5.0)
	#c = 1.0

	light = ("LIGHT_SPEED =" + str(c))
	masse = ("ME =" + str(me))
	massi = ("MI =" + str(mi))

	#files = ("runs = 'R_" + str(v) + "_c_" + str(c) + "'")

	folder = ("recon_008_MR_" + str(v) + "_c0_" + str(c))


	## create Warpx .pin file
	out = open('ssrecon_wv.pin','w')
	## write data into .pin file
	out.write('# -*- python -*- \n')
	out.write('# The following parameters has been randomly generated. \n')
	out.write(light)
	out.write("\n")
	out.write(masse)
	out.write("\n")
	out.write(massi)
	out.write("\n")
	#out.write(files)
	#out.write("\n")

	# add it to .pin file
	out.write('# -- End of randomly generated data. --')
	out.write(info)
	out.close

	os.mkdir(folder)
	os.system("mv ssrecon_wv.pin " + folder)
	os.system("cp cray.qsub " + folder)
