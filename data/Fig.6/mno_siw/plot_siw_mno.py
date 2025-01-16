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
w_norg_t2g, norg_t2g_re, norg_t2g_im = numpy.loadtxt("./norg/t2g/t2g.data", unpack = True, usecols = (0,1,2))
w_norg_eg, norg_eg_re, norg_eg_im = numpy.loadtxt("./norg/eg/eg.data", unpack = True, usecols = (0,1,2))
w_ctqmc_t2g, ctqmc_t2g_re, ctqmc_t2g_im = numpy.loadtxt("./ctqmc/t2g/t2g.data", unpack = True, usecols = (0,1,2))
w_ctqmc_eg, ctqmc_eg_re, ctqmc_eg_im = numpy.loadtxt("./ctqmc/eg/eg.data", unpack = True, usecols = (0,1,2))

plt.figure(0)
############
# t2g im part

# plot it
lines = plt.plot(w_norg_t2g, norg_t2g_im, w_ctqmc_t2g, ctqmc_t2g_im, alpha = 0.8, clip_on = True)

# setup line properties
plt.setp(lines[0], linewidth = 2.0, label = 'NORG')
plt.setp(lines[1], linewidth = 2.0, marker = 'o', label = 'CT-HYB')

# setup tics
plt.xticks()
plt.yticks([0.0,-0.5,-1.0,-1.5,-2.0,-2.5,-3.0,-3.5], ['0.0','-0.5','-1.0','-1.5','-2.0','-2.5','-3.0','-3.5'])
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
plt.ylim(-3.5,0.0)

# output the figure
plt.savefig("mno-mat-t2g-im.pdf",bbox_inches='tight')

plt.figure(1)
############
# eg im part

# plot it
lines = plt.plot(w_norg_eg, norg_eg_im, w_ctqmc_eg, ctqmc_eg_im, alpha = 0.8, clip_on = True)

# setup line properties
plt.setp(lines[0], linewidth = 2.0, label = 'NORG')
plt.setp(lines[1], linewidth = 2.0, marker = 'o', label = 'CT-HYB')

# setup tics
plt.xticks()
plt.yticks([0.0,-1.0,-2.0,-3.0,-4.0,-5.0], ['0.0','-1.0','-2.0','-3.0','-4.0','-5.0'])
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
plt.ylim(-5.0,0.0)

# output the figure
plt.savefig("mno-mat-eg-im.pdf",bbox_inches='tight')
