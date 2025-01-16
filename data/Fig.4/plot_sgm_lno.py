#!/usr/bin/env python

import sys
import numpy
import scipy.interpolate
import matplotlib
matplotlib.use("pdf")  # 设置后端
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon
from matplotlib.ticker import AutoMinorLocator
import matplotlib.gridspec as gridspec  # 导入gridspec模块
import datetime

# 设置字体和数学文本渲染
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['axes.linewidth'] = 1.00
plt.rcParams['lines.linewidth'] = 1.0
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}\usepackage{siunitx}\usepackage{newtxtext}\usepackage{newtxmath}'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def present():
    return datetime.datetime.now().strftime("%Y%m%d")

# 创建自定义网格布局
x_size = 8
fig = plt.figure(figsize=(x_size*1.618033988749895, x_size/2))
widths = [1, 1, 1, 1]  # 4列
gs = gridspec.GridSpec(2, 4, width_ratios=widths, wspace=0.26, hspace=0.0)

# 创建子图列表
axs = []
for i in range(4):
    axs.append(fig.add_subplot(gs[0, i]))  # 上排
for i in range(4):
    axs.append(fig.add_subplot(gs[1, i]))  # 下排

# left_annotate = (0.65, 0.82)
right_annotate = (0.03, 0.85)


a2d_annotate = (0.69, 0.90)
e2h_annotate = (0.75, 0.08)


# 读取Giw数据
w_norg_t2g_giw, norg_t2g_re_giw, norg_t2g_im_giw    = numpy.loadtxt("./lno_giw/norg-lno-dx.data", unpack=True, usecols=(0,1,2))
w_norg_eg_giw, norg_eg_re_giw, norg_eg_im_giw       = numpy.loadtxt("./lno_giw/norg-lno-dz.data", unpack=True, usecols=(0,1,2))
w_ctqmc_t2g_giw, ctqmc_t2g_re_giw, ctqmc_t2g_im_giw = numpy.loadtxt("./lno_giw/ctqmc-lno-dx.data", unpack=True, usecols=(0,1,2))
w_ctqmc_eg_giw, ctqmc_eg_re_giw, ctqmc_eg_im_giw    = numpy.loadtxt("./lno_giw/ctqmc-lno-dz.data", unpack=True, usecols=(0,1,2))

# 子图 (a)
axs[0].annotate(r'(a) $d_{x^2-y^2}$', a2d_annotate, xycoords='axes fraction', fontsize=12)
lines = axs[0].plot(w_norg_t2g_giw, norg_t2g_im_giw, w_ctqmc_t2g_giw, ctqmc_t2g_im_giw, alpha=0.8, clip_on=True)
plt.setp(lines[0], linewidth=2.0, label='NORG')
plt.setp(lines[1], linewidth=2.0, marker='o', ms=4, label='CT-HYB')
axs[0].set_ylabel(r"Im$G(i\omega_n)$", labelpad=3)
axs[0].set_xlim(0, 3)
axs[0].set_xticklabels([])
axs[0].set_ylim(-0.9, 0.0)
axs[0].set_yticks([ -0.6, -0.3, 0.0])
axs[0].tick_params(length=4, width=0.5, direction='in', top=False, right=False)


# 子图 (e)
axs[4].annotate(r'(e) $d_{z^2}$', e2h_annotate, xycoords='axes fraction', fontsize=12)
lines = axs[4].plot(w_norg_eg_giw, norg_eg_im_giw, w_ctqmc_eg_giw, ctqmc_eg_im_giw, alpha=0.8, clip_on=True)
plt.setp(lines[0], linewidth=2.0, label='NORG')
plt.setp(lines[1], linewidth=2.0, marker='o', ms=4, label='CT-HYB')
axs[4].set_xlabel(r"$\omega_n$", labelpad=-1.0)
axs[4].set_ylabel(r"Im$G(i\omega_n)$", labelpad=3)
# axs[0].legend(frameon=False, loc='upper left', ncol=1)
axs[4].set_xlim(0, 3)
axs[4].set_xticks([0, 1, 2, 3])
axs[4].set_ylim(-1.5, 0.0)
axs[4].set_yticks([-1.5, -1.0, -0.5, 0.0])
axs[4].tick_params(length=4, width=0.5, direction='in', top=True, right=False)


# 读取Siw数据
w_norg_t2g_siw, norg_t2g_re_siw, norg_t2g_im_siw    = numpy.loadtxt("./lno_siw/norg/dx2-y2/dx2-y2.data", unpack=True, usecols=(0,1,2))
w_norg_eg_siw, norg_eg_re_siw, norg_eg_im_siw       = numpy.loadtxt("./lno_siw/norg/dz2/dz2.data", unpack=True, usecols=(0,1,2))
w_ctqmc_t2g_siw, ctqmc_t2g_re_siw, ctqmc_t2g_im_siw = numpy.loadtxt("./lno_siw/ctqmc/dx2-y2/dx2-y2.data", unpack=True, usecols=(0,1,2))
w_ctqmc_eg_siw, ctqmc_eg_re_siw, ctqmc_eg_im_siw    = numpy.loadtxt("./lno_siw/ctqmc/dz2/dz2.data", unpack=True, usecols=(0,1,2))

# 子图 (b)
axs[1].annotate(r'(b) $d_{x^2-y^2}$', a2d_annotate, xycoords='axes fraction', fontsize=12)
lines = axs[1].plot(w_norg_t2g_siw, norg_t2g_im_siw, w_ctqmc_t2g_siw, ctqmc_t2g_im_siw, alpha=0.8, clip_on=True)
plt.setp(lines[0], linewidth=2.0, label='NORG')
plt.setp(lines[1], linewidth=2.0, marker='o', ms=4, label='CT-HYB')
# axs[1].set_xlabel(r"$\omega_n$")
axs[1].set_ylabel(r"Im$\Sigma(i\omega_n)$", labelpad=3)
axs[1].set_xlim(0, 8)
axs[1].set_xticklabels([])
axs[1].set_ylim(-0.8, 0.0)
axs[1].set_xticks([0, 4, 8])
axs[1].set_yticks([  -0.6, -0.3, 0.0])
axs[1].tick_params(length=4, width=0.5, direction='in', top=False, right=True)


# 子图 (f)
axs[5].annotate(r'(f) $d_{z^2}$', e2h_annotate, xycoords='axes fraction', fontsize=12)
lines = axs[5].plot(w_norg_eg_siw, norg_eg_im_siw, w_ctqmc_eg_siw, ctqmc_eg_im_siw, alpha=0.8, clip_on=True)
plt.setp(lines[0], linewidth=2.0, label='NORG')
plt.setp(lines[1], linewidth=2.0, marker='o', ms=4, label='CT-HYB')
axs[5].set_xlabel(r"$\omega_n$", labelpad=-1.0)
axs[5].set_ylabel(r"Im$\Sigma(i\omega_n)$", labelpad=3)
axs[5].set_xlim(0, 8)
axs[5].set_xticks([0, 4, 8])
axs[5].set_ylim(-1.0, 0.0)
axs[5].set_yticks([-1.0, -0.5, 0])
axs[5].tick_params(length=4, width=0.5, direction='in', top=True, right=True)

# ----------------------------------------------------------------------------------------------------------------------------------
# 读取Sigma数据
w_norg_t2g_sigma, norg_t2g_re_sigma, norg_t2g_im_sigma      = numpy.loadtxt("./lno-sigma/norg-dx2-sigma.data", unpack = True, usecols = (0,1,2))
w_norg_eg_sigma, norg_eg_re_sigma, norg_eg_im_sigma         = numpy.loadtxt("./lno-sigma/norg-dz2-sigma.data", unpack = True, usecols = (0,1,2))
w_ctqmc_t2g_sigma, ctqmc_t2g_re_sigma, ctqmc_t2g_im_sigma   = numpy.loadtxt("./lno-sigma/ctqmc-dx2-sigma.data", unpack = True, usecols = (0,1,2))
w_ctqmc_eg_sigma, ctqmc_eg_re_sigma, ctqmc_eg_im_sigma      = numpy.loadtxt("./lno-sigma/ctqmc-dz2-sigma.data", unpack = True, usecols = (0,1,2))
cdgh_annotate = (0.05, 0.1)

# 子图 (c)
axs[2].annotate(r'(c) $d_{x^2-y^2}$', a2d_annotate, xycoords='axes fraction', fontsize=12)
lines = axs[2].plot(w_norg_t2g_sigma, norg_t2g_re_sigma, w_ctqmc_t2g_sigma, ctqmc_t2g_re_sigma, alpha=0.8, clip_on=True)
axs[2].plot([0, 0], [0, 7], 'k--')
plt.setp(lines[0], linewidth=2.0, label='NORG')
plt.setp(lines[1], linewidth=2.0, label='CT-HYB')
axs[2].set_ylabel(r"Re$\Sigma(\omega)$", labelpad=3)
axs[2].set_xlim(-6, 4)
axs[2].set_xticks([-6, -4, -2, 0, 2, 4])
axs[2].set_xticklabels([])
axs[2].set_ylim(4, 7)
axs[2].set_yticks([ 5, 6, 7])
axs[2].tick_params(length=4, width=0.5, direction='in', top=False, right=False)


# 子图 (g)
axs[6].annotate(r'(g) $d_{z^2}$', e2h_annotate, xycoords='axes fraction', fontsize=12)
lines = axs[6].plot(w_norg_eg_sigma, norg_eg_re_sigma, w_ctqmc_eg_sigma, ctqmc_eg_re_sigma, alpha=0.8, clip_on=True)
axs[6].plot([0, 0], [0, 8], 'k--')
plt.setp(lines[0], linewidth=2.0, label='NORG')
plt.setp(lines[1], linewidth=2.0, label='CT-HYB')
axs[6].set_xlabel(r"$\omega$ (eV)", labelpad= 0.0)
axs[6].set_ylabel(r"Re$\Sigma(\omega)$", labelpad=3)
# axs[2].legend(frameon=False, loc='upper left')
axs[6].set_xlim(-6, 4)
axs[6].set_xticks([-6, -4, -2, 0, 2, 4])
axs[6].set_ylim(4.0,7.0)
axs[6].set_yticks([4, 5, 6, 7])
axs[6].tick_params(length=4, width=0.5, direction='in', top=True, right=False)


# 子图 (d)
axs[3].annotate(r'(d) $d_{x^2-y^2}$', a2d_annotate, xycoords='axes fraction', fontsize=12)
lines = axs[3].plot(w_norg_t2g_sigma, norg_t2g_im_sigma, w_ctqmc_t2g_sigma, ctqmc_t2g_im_sigma, alpha=0.8, clip_on=True)
axs[3].plot([0, 0], [0, -4], 'k--')
plt.setp(lines[0], linewidth=2.0)
plt.setp(lines[1], linewidth=2.0)
# axs[3].set_xlabel(r"$\omega$ (eV)")
axs[3].set_ylabel(r"Im$\Sigma(\omega)$", labelpad=3)
axs[3].set_xlim(-6, 4)
# axs[5].set_xticks([-4, -2, 0, 2, 4, 6, 8])
axs[3].set_ylim(-2.5, 0.0)
axs[3].set_yticks([-2, -1, 0])
axs[3].set_xticklabels([])
axs[3].tick_params(length=4, width=0.5, direction='in', top=False, right=True)

# 子图 (h)
axs[7].annotate(r'(h) $d_{z^2}$', e2h_annotate, xycoords='axes fraction', fontsize=12)
lines = axs[7].plot(w_norg_eg_sigma, norg_eg_im_sigma, w_ctqmc_eg_sigma, ctqmc_eg_im_sigma, alpha=0.8, clip_on=True)
axs[7].plot([0, 0], [0, -4], 'k--')
plt.setp(lines[0], linewidth=2.0)
plt.setp(lines[1], linewidth=2.0)
axs[7].set_xlabel(r"$\omega$ (eV)", labelpad= 0.0)
axs[7].set_ylabel(r"Im$\Sigma(\omega)$", labelpad=3)
axs[7].set_xlim(-6, 4)
# axs[7].set_xticks([-4, -2, 0, 2, 4, 6, 8])
axs[7].set_ylim(-4.0, 0.0)
axs[7].set_yticks([-4, -2, 0])
axs[7].tick_params(length=4, width=0.5, direction='in', top=True, right=True)


axs[0].legend(frameon=False, loc='lower left', ncol=1, bbox_to_anchor=(-0.02, -0.04))
axs[2].legend(frameon=False, loc='lower left', ncol=1, bbox_to_anchor=(-0.02, -0.04))

# 添加次要刻度
for i, ax in enumerate(axs):
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    if i >= 4:  # 针对子图 e、f、g、h
        ax.tick_params(which='minor', length=2, width=0.5, top=True)
    else:
        ax.tick_params(which='minor', length=2, width=0.5)
# 调整布局并保存图片
plt.tight_layout()
plt.savefig("lno-sgm"+present(), bbox_inches='tight', pad_inches=0.00625, transparent=True)
plt.savefig("../../tex-zen/lno-sgm", bbox_inches='tight', pad_inches=0.00625, transparent=True)