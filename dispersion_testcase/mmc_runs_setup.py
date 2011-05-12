import os
import random
import string
	

# get paramenters that remain the same for all runs
f = open('input3.py', 'r')
info = f.read()
f.close

#write all the values on a file
avg = open('mmc_values.txt','w')

# Number of runs
last = 36  # number of runs

for count in range(last):
	# uniform varying values of mass ration from
	# a to b given a certain value of ME
	v = random.uniform(1.0e-10,1.0e-4)
	#v = random.gauss(3,1.0)

	ux = ("value = " + str(v))

	folder = ("advect_001_U_" + str(v))


	## create Warpx .pin file
	out = open('advect_0.pin','w')
	## write data into .pin file
	out.write('# -*- python -*- \n')
	out.write('# The following parameters has been randomly generated. \n')
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
	out.close

	os.mkdir(folder)
	os.system("mv advect_0.pin " + folder)
	avg.write(str(v))
	avg.write('\n')

	if count <= round(last/2):
		out = open('advect_1.pin','w')
		out.write('# -*- python -*- \n')
		out.write('# The following parameters has been randomly generated. \n')
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
		out.close
		os.system("mv advect_1.pin " + folder)

	if count <= round(last/4):
		out = open('advect_2.pin','w')
		out.write('# -*- python -*- \n')
		out.write('# The following parameters has been randomly generated. \n')
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
		out.close
		os.system("mv advect_2.pin " + folder)

avg.close
