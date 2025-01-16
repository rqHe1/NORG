#!/usr/bin/env python

import sys
import numpy
import scipy.interpolate
import matplotlib
# matplotlib.use("Agg") # setup backend
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

# Adjusted figure size
x_size = 2.5  
fig, ax = plt.subplots(1, 1, figsize=(x_size*1.618033988749895, x_size))  
plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)

colors = {
    "NORG_t2g": "tab:blue",
    "NORG_eg": "tab:purple",
    "CTQMC_t2g": "tab:green",
    "CTQMC_eg": "tab:red",
    "Exp": "tab:orange"
}

# Read data
w_norg_t2g, norg_t2g = numpy.loadtxt("./data-etaE-1/norg-t2g-Aout.data", unpack=True, usecols=(0, 1))
w_norg_eg, norg_eg = numpy.loadtxt("./data-etaE-1/norg-eg-Aout.data", unpack=True, usecols=(0, 1))
w_ctqmc_t2g, ctqmc_t2g = numpy.loadtxt("./data-etaE-1/ctqmc-t2g-Aout.data", unpack=True, usecols=(0, 1))
w_ctqmc_eg, ctqmc_eg = numpy.loadtxt("./data-etaE-1/ctqmc-eg-Aout.data", unpack=True, usecols=(0, 1))

w_exp, y_exp = numpy.loadtxt("./data-etaE-1/svo-exp.data", unpack=True, usecols=(0, 1))
y_exp = y_exp * 0.01
# Plot data with specified colors
lines1, = ax.plot(
    w_norg_t2g, norg_t2g,
    color=colors["NORG_t2g"],
    alpha=0.8,
    label=r"$t_{2g}$ NORG"
)
lines3, = ax.plot(
    w_ctqmc_t2g, ctqmc_t2g,
    color=colors["CTQMC_t2g"],
    alpha=0.8,
    label=r"$t_{2g}$ CT-HYB"
)
lines2, = ax.plot(
    w_norg_eg, norg_eg,
    color=colors["NORG_eg"],
    alpha=0.8,
    label=r"$e_{g}$ NORG"
)
lines4, = ax.plot(
    w_ctqmc_eg, ctqmc_eg,
    color=colors["CTQMC_eg"],
    alpha=0.8,
    label=r"$e_{g}$ CT-HYB"
)

# EXP
exp_line, = ax.plot(
    w_exp, y_exp,
    color=colors["Exp"],        
    marker='o',                 
    linestyle='-',              
    markersize=3,               
    alpha=0.7,                  
    label='Expt.',              
    markerfacecolor='none',     
    zorder=5                    
)


# Add vertical line at omega=0
ax.axvline(0, linestyle='--', color='black')

# Setup ticks
ax.set_yticks([0, 0.4, 0.8, 1.2])
ax.tick_params(length=4, width=0.5, direction='in', top=False, right=False)
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.tick_params(which='minor', length=2, width=0.5, direction='in', top=False, right=False)

# Setup labels
ax.set_ylabel(r'$A(\omega)$')
ax.set_xlabel(r'$\omega$ (eV)', labelpad=0.5)

# Add legend
ax.legend(frameon=False, loc='upper right')

# Setup x and y range
ax.set_xlim(-4, 8)
ax.set_ylim(0, 1.1)

# for i, ax in enumerate(axs):
#     ax.xaxis.set_minor_locator(AutoMinorLocator(2))
#     ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#     if i >= 4:  # 针对子图 e、f、g、h
#         ax.tick_params(which='minor', length=2, width=0.5, top=True)
#     else:
#         ax.tick_params(which='minor', length=2, width=0.5)
# Add overall title or annotation if needed
# ax.annotate('Combined Plot of NORG and CT-HYB', xy=(0.5, 1.05), xycoords='axes fraction', ha='center', fontsize=12)

# Output the figure
plt.savefig("svo-dos"+present()+".pdf", bbox_inches='tight', pad_inches=0.00625, transparent=True)
plt.savefig("../../tex-zen/svo-dos.pdf", bbox_inches='tight', pad_inches=0.00625, transparent=True)
plt.close()