# Reconnected magnetic flux for one run folder.
#
# Run inside a run folder. For each WARPX output prefix, reads frames
# 0..frames of 'qnew', integrates |B_y| (component 14) along the mid-plane
# y = Ly/2, normalises the time series so that flux(0) = 2*wci (to match the
# GEM challenge), and writes it to <prefix>.dat (one value per frame).
#
# Usage:
#   python calc_flux.py ssrecon_wv_0 ...  # the prefixes given
#   python calc_flux.py                   # every prefix found in the folder:
#                                         # ssrecon_wv_0/1/2 (MC/MMC levels)
#                                         # and recon_pcm (PCM)
#
# Replaces calc_flux_0.py, calc_flux_1.py, calc_flux_2.py, calc_flux_mc.py,
# calc_flux_mmc.py and calc_flux_pcm.py, which did the same computation
# for fixed prefixes.
import os
import sys
import numpy
import wxdata as wxdata2

frames = 40
Ly = 12.8
wci = 0.1

KNOWN = ['ssrecon_wv_0', 'ssrecon_wv_1', 'ssrecon_wv_2', 'recon_pcm']


def reconnected_flux(prefix):
    flux = numpy.zeros(frames+1, float)
    for n in range(0, frames+1):
        dh = wxdata2.WxData(prefix, n)
        q = dh.read('qnew')
        dx = (q.grid.upperBounds[1]-q.grid.lowerBounds[1])/q.grid.numPhysCells[1]
        ny = q.grid.numPhysCells[1]//2
        by = q[:, ny, 14]
        flux[n] = dx*numpy.sum(numpy.fabs(by))
        dh.close()
    flux = flux/(2*Ly)
    flux = 2*wci*flux/flux[0]  # rescale to match GEM conditions
    return flux


if __name__ == '__main__':
    prefixes = sys.argv[1:]
    if not prefixes:
        prefixes = [p for p in KNOWN if os.path.exists(p + '.pin')]
    if not prefixes:
        sys.exit('calc_flux.py: no prefix given and no known .pin file found')
    for prefix in prefixes:
        numpy.savetxt(prefix + '.dat', reconnected_flux(prefix))
