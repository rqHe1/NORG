import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.use("pdf") # setup backend
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import AutoMinorLocator

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




# Read data
w_norg_t2g_1, norg_t2g_1      = np.loadtxt("./mno-1.00/norg-dos/t2g/Aout.data", unpack=True, usecols=(0,1))
w_norg_eg_1, norg_eg_1        = np.loadtxt("./mno-1.00/norg-dos/eg/Aout.data", unpack=True, usecols=(0,1))
w_ctqmc_t2g_1, ctqmc_t2g_1    = np.loadtxt("./mno-1.00/ctqmc-dos/t2g/Aout.data", unpack=True, usecols=(0,1))
w_ctqmc_eg_1, ctqmc_eg_1      = np.loadtxt("./mno-1.00/ctqmc-dos/eg/Aout.data", unpack=True, usecols=(0,1))

w_norg_t2g_1    = w_norg_t2g_1  + 1.7
w_norg_eg_1     = w_norg_eg_1   + 1.7
w_ctqmc_t2g_1   = w_ctqmc_t2g_1 + 1.7
w_ctqmc_eg_1    = w_ctqmc_eg_1  + 1.7

w_norg_t2g_2, norg_t2g_2      = np.loadtxt("./mno-0.53/norg-t2g-Aout.data", unpack=True, usecols=(0,1))
w_norg_eg_2, norg_eg_2        = np.loadtxt("./mno-0.53/norg-eg-Aout.data", unpack=True, usecols=(0,1))
w_ctqmc_t2g_2, ctqmc_t2g_2    = np.loadtxt("./mno-0.53/ctqmc-t2g-Aout.data", unpack=True, usecols=(0,1))
w_ctqmc_eg_2, ctqmc_eg_2      = np.loadtxt("./mno-0.53/ctqmc-eg-Aout.data", unpack=True, usecols=(0,1))



w_exp_1, y_exp_1 = np.loadtxt("./mno-1.00/mno-exp-lower.data", unpack=True, usecols=(0, 1))
w_exp_2, y_exp_2 = np.loadtxt("./mno-1.00/mno-exp-upper.data", unpack=True, usecols=(0, 1))
y_exp_1 = y_exp_1 * 0.001
y_exp_2 = y_exp_2 * 0.001
w_exp_1 = w_exp_1 + 0.2
w_exp_2 = w_exp_2 + 0.2
# Create the main figure with two subplots
fig = plt.figure(figsize=(4.0, 5))
gs_lower = GridSpec(2, 2, wspace=0, hspace=0, top=0.55, bottom=0.1 )



# Plot for Two-orbital Wide band case
ax1 = fig.add_subplot(gs_lower[0,0])
ax1.plot(w_norg_t2g_1, norg_t2g_1, alpha = 0.8, clip_on = True)#, label = r"$t_{2g}$")
ax1.plot(w_norg_eg_1, norg_eg_1  , alpha = 0.8, clip_on = True)#, label = r"$e_{g}$" )
ax1.plot( w_exp_1[::4], y_exp_1[::4], marker='o',linestyle='-',markersize=3,alpha=0.7,label='XPS',markerfacecolor='none',zorder=5)
ax1.plot( w_exp_2[::4], y_exp_2[::4], marker='o',linestyle='-',markersize=3,alpha=0.7,label='BIS',markerfacecolor='none',zorder=5)
ax1.set_xlim(-8, 6)
# ax1.set_xticks([-0.2, 0.0, 0.2, 0.4])
ax1.set_ylim(0, 0.6)
ax1.set_yticks([0.0, 0.2, 0.4, 0.6])
ax1.set_xticks([-8, -4, 0, 4])
ax1.set_xticklabels([])
ax1.set_ylabel(r'$A(\omega)$', labelpad=3)
ax1.text(0.02, 0.88, '(a) NORG', transform=ax1.transAxes)
# ax1.legend(frameon=False, loc='upper right', bbox_to_anchor=(0.63, 0.9), prop={'size': 8})#, handlelength=1.0)
# ax1.legend(frameon=False, loc='upper right', bbox_to_anchor=(1, 1.06), prop={'size': 8})#, handlelength=1.0)
ax1.legend(frameon=False, loc='upper left', bbox_to_anchor=(0, 0.9)) 
# ax1.legend(frameon=False, loc='upper left', bbox_to_anchor=(-0.02, -0.15), prop={'size': 8})#, handlelength=1.0)



ax2 = fig.add_subplot(gs_lower[1,0])
ax2.plot(w_ctqmc_t2g_1, ctqmc_t2g_1, alpha = 0.8, clip_on = True, label = r"$t_{2g}$")
ax2.plot(w_ctqmc_eg_1, ctqmc_eg_1,   alpha = 0.8, clip_on = True, label = r"$e_{g}$" )
ax2.plot( w_exp_1[::4], y_exp_1[::4], marker='o',linestyle='-',markersize=3,alpha=0.7,label='XPS',markerfacecolor='none',zorder=5)
ax2.plot( w_exp_2[::4], y_exp_2[::4], marker='o',linestyle='-',markersize=3,alpha=0.7,label='BIS',markerfacecolor='none',zorder=5)
ax2.set_xlim(-8, 6)
ax2.set_ylim(0, 0.6)
ax2.set_yticks([0.0, 0.2, 0.4])
ax2.set_xticks([-8, -4, 0, 4])
ax2.text(0.02, 0.88, '(c) CT-HYB', transform=ax2.transAxes)
ax2.set_xlabel(r'$\omega$', labelpad=0.20)
ax2.set_ylabel(r'$A(\omega)$', labelpad=3)
# ax1.legend(frameon=False, loc='upper right', bbox_to_anchor=(0.63, 0.9), prop={'size': 8})#, handlelength=1.0)
# ax1.legend(frameon=False, loc='upper right', bbox_to_anchor=(1, 1.06), prop={'size': 8})#, handlelength=1.0)





# Plot for Two-orbital Narrow band case
ax3 = fig.add_subplot(gs_lower[0,1])
ax3.plot(w_norg_t2g_2, norg_t2g_2, label = r"$t_{2g}$", alpha = 0.8, clip_on = True)
ax3.plot(w_norg_eg_2, norg_eg_2  , label = r"$e_{g}$" ,   alpha = 0.8, clip_on = True)
ax3.set_xlim(-8, 6)
ax3.set_ylim(0, 0.6)
# ax3.set_xticks([-0.1, 0.1, 0.3, 0.5])
ax3.set_yticks([0.0, 0.2, 0.4])
ax3.set_xticks([-8, -4, 0, 4])
ax3.set_yticklabels([])
ax3.set_xticklabels([])
ax3.text(0.02, 0.88, '(b) NORG', transform=ax3.transAxes)
# ax3.legend(frameon=False, loc='upper right', bbox_to_anchor=(0.63, 0.9), prop={'size': 8})#, handlelength=1.0)
# ax3.legend(frameon=False, loc='upper left', bbox_to_anchor=(0.0, 0.8), prop={'size': 8})#, handlelength=1.0)
# ax3.legend(frameon=False, loc='upper right', bbox_to_anchor=(1, 1), prop={'size': 8}) 
ax3.legend(frameon=False, loc='upper left', bbox_to_anchor=(0, 0.9)) 

ax4 = fig.add_subplot(gs_lower[1,1])
ax4.plot(w_ctqmc_t2g_2, ctqmc_t2g_2, alpha = 0.8, clip_on = True)
ax4.plot(w_ctqmc_eg_2, ctqmc_eg_2,   alpha = 0.8, clip_on = True)
ax4.set_xlim(-8, 6)
ax4.set_ylim(0, 0.6)
ax4.set_yticklabels([])
# ax4.set_ylim([0, 0.25])
ax4.set_yticks([0.0, 0.2, 0.4])
ax4.set_xticks([-8, -4, 0, 4])
ax4.text(0.02, 0.88, '(d) CT-HYB', transform=ax4.transAxes)
ax4.set_xlabel(r'$\omega$', labelpad=0.20)
# ax4.legend(frameon=False, loc='upper left', bbox_to_anchor=(-0.02, 0.9), prop={'size': 8})#, handlelength=1.0)

# ax3.text(0.5, 1.03, 'MnO-53%% cell',ha='center', transform=ax3.transAxes)
ax1.text(0.5, 1.05, 'MnO',ha='center', transform=ax1.transAxes)
ax3.text(0.5, 1.05, "MnO compressed", ha='center', transform=ax3.transAxes)
plt.subplots_adjust(top=1.1)
# ax.text(0.5, 0.6, 'MnO-53%% cell', ha='center')

axes = [ax1, ax2, ax3, ax4]
for ax in axes:
    ax.axvline(0, linestyle=(0, (2, 2)), color='black', linewidth=0.8)
    ax.tick_params(length=4, width=0.5, direction='in', top=False, right=False)
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which='minor', length=2, width=0.5, direction='in', top=False, right=False)

# Save the figure
plt.tight_layout()
plt.savefig("mno-dos"+present()+".pdf", format='pdf', bbox_inches='tight', pad_inches=0.00625, transparent=True)
plt.savefig("../../tex-zen/mno-dos.pdf", format='pdf', bbox_inches='tight', pad_inches=0.00625, transparent=True)
