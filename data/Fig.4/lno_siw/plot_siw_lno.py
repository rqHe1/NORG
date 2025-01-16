#!/usr/bin/env python

import sys
import numpy
import scipy.interpolate
import matplotlib
#matplotlib.use("Agg") # setup backend
matplotlib.use("pdf") # setup backend
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon

matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'

# read data
w_norg_dx2, norg_dx2_re, norg_dx2_im = numpy.loadtxt("./norg/dx2-y2/dx2-y2.data", unpack = True, usecols = (0,1,2))
w_norg_dz2, norg_dz2_re, norg_dz2_im = numpy.loadtxt("./norg/dz2/dz2.data", unpack = True, usecols = (0,1,2))
w_ctqmc_dx2, ctqmc_dx2_re, ctqmc_dx2_im = numpy.loadtxt("./ctqmc/dx2-y2/dx2-y2.data", unpack = True, usecols = (0,1,2))
w_ctqmc_dz2, ctqmc_dz2_re, ctqmc_dz2_im = numpy.loadtxt("./ctqmc/dz2/dz2.data", unpack = True, usecols = (0,1,2))

plt.figure(0)
############
# dx2 im part

# plot it
lines = plt.plot(w_norg_dx2, norg_dx2_im, w_ctqmc_dx2, ctqmc_dx2_im, alpha = 0.8, clip_on = True)

# setup line properties
plt.setp(lines[0], linewidth = 2.0, label = 'NORG')
plt.setp(lines[1], linewidth = 2.0, marker = 'o', label = 'CT-HYB')

# setup tics
plt.xticks()
plt.yticks([0.0,-0.2,-0.4,-0.6,-0.8],['0.00', '-0.20', '-0.40', '-0.60', '-0.80'])
plt.tick_params(length = 8, width = 1.0, which = 'major', direction = 'in')
plt.tick_params(length = 4, width = 0.5, which = 'minor', direction = 'in')
plt.gca().xaxis.set_ticks_position('both')
plt.gca().yaxis.set_ticks_position('both')

# setup labels
plt.xlabel(r"$i\omega_n$")
plt.ylabel(r"Im$\Sigma(i\omega_n)$")
plt.legend(framealpha = 0.5, loc = 'lower right', ncol = 1)

# setup x and y range
plt.xlim(0,8)
plt.ylim(-0.8,0.0)

# output the figure
plt.savefig("lno-mat-dx2-y2-im.pdf",bbox_inches='tight')

plt.figure(1)
############
# eg im part

# plot it
lines = plt.plot(w_norg_dz2, norg_dz2_im, w_ctqmc_dz2, ctqmc_dz2_im, alpha = 0.8, clip_on = True)

# setup line properties
plt.setp(lines[0], linewidth = 2.0, label = 'NORG')
plt.setp(lines[1], linewidth = 2.0, marker = 'o', label = 'CT-HYB')

# setup tics
plt.xticks()
plt.yticks([0.0,-0.2,-0.4,-0.6,-0.8,-1.0],['0.00', '-0.20', '-0.40', '-0.60', '-0.80', '-1.00'])
plt.tick_params(length = 8, width = 1.0, which = 'major', direction = 'in')
plt.tick_params(length = 4, width = 0.5, which = 'minor', direction = 'in')
plt.gca().xaxis.set_ticks_position('both')
plt.gca().yaxis.set_ticks_position('both')

# setup labels
plt.xlabel(r"$i\omega_n$")
plt.ylabel(r"Im$\Sigma(i\omega_n)$")
plt.legend(framealpha = 0.5, loc = 'lower right', ncol = 1)

# setup x and y range
plt.xlim(0,8)
plt.ylim(-1,0.0)

# output the figure
plt.savefig("lno-mat-dz2-im.pdf",bbox_inches='tight')
