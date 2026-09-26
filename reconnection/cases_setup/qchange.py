# Cancel a range of batch jobs with qdel.
#
# Usage: python qchange.py FIRST LAST   (cancels job IDs FIRST .. LAST-1)
import os
import sys

first, last = int(sys.argv[1]), int(sys.argv[2])

for k in range(first, last):
    print(k)
    os.system("qdel "+str(k))
