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
from matplotlib.ticker import AutoMinorLocator

matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'

import datetime
def present():
    return datetime.datetime.now().strftime("%Y%m%d")
# Set general font properties
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 7  # Adjusted font size for better compatibility with article text
plt.rcParams['axes.linewidth'] = 1.00  # Adjusted axis line width for better compatibility with article text
plt.rcParams['lines.linewidth'] = 1.0  # Sets default data line width
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}\usepackage{siunitx}\usepackage{newtxtext}\usepackage{newtxmath}'

# Set font family for all text elements
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# Adjusted figure size in subplots
x_size = 2.5
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(x_size*1.618033988749895, x_size))  # Adjusted size for smaller output
plt.subplots_adjust(hspace=0.0)

# Read data
w_norg_t2g, norg_t2g = numpy.loadtxt("./data-etaE-1/norg-dx2-Aout.data", unpack = True, usecols = (0,1))
w_norg_eg, norg_eg = numpy.loadtxt("./data-etaE-1/norg-dz2-Aout.data", unpack = True, usecols = (0,1))
w_ctqmc_t2g, ctqmc_t2g = numpy.loadtxt("./data-etaE-1/ctqmc-dx2-Aout.data", unpack = True, usecols = (0,1))
w_ctqmc_eg, ctqmc_eg = numpy.loadtxt("./data-etaE-1/ctqmc-dz2-Aout.data", unpack = True, usecols = (0,1))

# plot it
# fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
# plt.subplots_adjust(hspace = 0.0)
lines1 = ax1.plot(w_norg_t2g, norg_t2g, alpha = 0.8, clip_on = True, label = 'NORG')
lines3 = ax1.plot(w_ctqmc_t2g, ctqmc_t2g, alpha = 0.8, clip_on = True, label = 'CT-HYB')
lines2 = ax2.plot(w_norg_eg, norg_eg, alpha = 0.8, clip_on = True, label = 'NORG')
lines4 = ax2.plot(w_ctqmc_eg, ctqmc_eg, alpha = 0.8, clip_on = True, label = 'CT-HYB')
# ax1.fill_between(w_norg_t2g, norg_t2g, color = lines1[0].get_color(), alpha = 0.1)
# ax1.fill_between(w_norg_eg, norg_eg, color = lines2[0].get_color(), alpha = 0.1)
# ax2.fill_between(w_ctqmc_t2g, ctqmc_t2g, color = lines3[0].get_color(), alpha = 0.1)
# ax2.fill_between(w_ctqmc_eg, ctqmc_eg, color = lines4[0].get_color(), alpha = 0.1)
ax1.axvline(0, linestyle = '--', color = 'black')
ax2.axvline(0, linestyle = '--', color = 'black')

# setup line properties

# setup tics
# ax1.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
# ax2.set_yticks([0, 0.2, 0.4, 0.6,])

ax1.set_ylim(0, 0.3)
ax1.set_yticks([0, 0.1, 0.2, 0.3])
ax2.set_ylim(0, 0.3)
ax2.set_yticks([0, 0.1, 0.2])
ax1.tick_params(length = 4, width = 0.5, direction = 'in', top=False, right=False)
ax2.tick_params(length = 4, width = 0.5, direction = 'in', top=True, right=False)

# setup labels
ax2.set_ylabel(r'$A(\omega)$')
ax2.set_xlabel(r'$\omega$ (eV)', labelpad=0.5)
ax1.set_ylabel(r'$A(\omega)$')
ax1.annotate(r'(a) $d_{x^2-y^2}$', xy=(0.02, 0.85),  xycoords='axes fraction')
ax2.annotate(r'(b) $d_{z^2}$', xy=(0.02, 0.85), xycoords='axes fraction')
# ax2.legend()
ax1.legend(frameon=False, loc='upper right')

# setup x and y range
ax1.set_xlim(-6, 4)
ax2.set_xlim(-6, 4)


ax1.tick_params(length=4, width=0.5, direction='in', top=False, right=False)
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
ax1.tick_params(which='minor', length=2, width=0.5, direction='in', top=False, right=False)

ax2.tick_params(length=4, width=0.5, direction='in', top=True, right=False)
ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
ax2.tick_params(which='minor', length=2, width=0.5, direction='in', top=True, right=False)




# Output the figure
plt.savefig("lno-dos"+present(), bbox_inches='tight', pad_inches=0.00625, transparent=True)
plt.savefig("../../tex-zen/lno-dos", bbox_inches='tight', pad_inches=0.00625, transparent=True)