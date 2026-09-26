# Convergence of plain Monte Carlo with the number of samples.
#
# Draws random subsets (with replacement) of increasing size from a large set
# of finished MC runs in ./MC_RUNS/advect_002_U_*, and for each subset size
# computes the sample mean and unbiased sample variance of component 1 of
# 'qnew' at the last output frame.
#
# Writes mc_stat_<N>.dat for each subset size N, with two columns: mean and
# variance at each grid cell.
import glob
import os
import numpy
import wxdata as wxdata2

frame = 10    # output frame to analyse
d = sorted(glob.glob("./MC_RUNS/advect_002_U_*"))

Num = [36,72,150,520,750,1150,1500,2500,3500,4500]

for N in Num:
    mn = numpy.random.randint(len(d), size=N)

    meanf = 0.
    varf  = 0.
    for m in mn:
        filename = os.path.join(d[m], "advect_mc")
        dh = wxdata2.WxData(filename, frame)
        q = dh.read('qnew')
        meanf = meanf + q[:,1]
        varf  = varf + q[:,1]*q[:,1]
        dh.close()

    varf  = varf/(N-1.) - meanf*meanf/(N*N-N)
    meanf = meanf/N
    numpy.savetxt('mc_stat_%d.dat' % N, numpy.column_stack((meanf, varf)))
