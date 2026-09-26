# Set the initial velocity in input3.inp by replacing the placeholder
# "u_i = value", e.g. to run a single deterministic case by hand.
#
# Usage: python replace_text.py [u_i]      (default u_i = 1.0e-5)
# Edits input3.inp in place; keep a copy if you still need the template.
import sys
import fileinput

file = "input3.inp"
u_i = sys.argv[1] if len(sys.argv) > 1 else "1.0e-5"

for line in fileinput.input(file, inplace=1):
    if "        u_i = value" in line:
        line = line.replace("        u_i = value","        u_i = "+u_i)
    sys.stdout.write(line)
