import sys
import fileinput

file = "input3.inp"

for line in fileinput.input(file, inplace=1):
    if "        u_i = value" in line:
        line = line.replace("        u_i = value","        u_i = 1.0e-5")
		sys.stdout.write(line)
