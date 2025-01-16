#!/usr/local/bin/python3
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

filename = 'svo'
threshold = 0.20  # 根据你的数据调整这个阈值

def fit_near_zero(x_data, y_data, threshold):
    near_zero_indices = np.abs(x_data) < threshold
    x_near_zero = x_data[near_zero_indices]
    y_near_zero = y_data[near_zero_indices]
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_near_zero, y_near_zero)
    return slope, intercept

def plot_fit(x_data, y_data, fit_x, fit_y, slope, label, color):
    plt.figure()  # 每次绘制一个新图
    plt.scatter(x_data, y_data, label=label)
    plt.plot(fit_x, fit_y, '--', label=f'Fit {label}: $m^*$: {1-slope:.4f}', color=color)
    plt.legend()
    plt.title(filename)
    plt.xlabel(r'$\omega$ (eV)')
    plt.ylabel(r'Re[$\Sigma(\omega)$]')
    plt.xlim(-3, 3)
    plt.ylim(min(y_data) - 0.1 * np.ptp(y_data), max(y_data) + 0.1 * np.ptp(y_data))
    plt.savefig(filename + f'_{label}_Refit.pdf', format='pdf')
    # plt.savefig(filename + f'_{label}_Refit.png')
    plt.close()  # 关闭图以便下次绘制新的图


# 分别提取x和y数据
fit_x = np.linspace(-10, 10, 1000)

x_norg_data, norg_eg_data = np.loadtxt("./norg-eg-sigma.data", unpack=True, usecols=(0,1))
slope, intercept = fit_near_zero(x_norg_data, norg_eg_data, threshold)
fit_y = slope * fit_x + intercept
plot_fit(x_norg_data, norg_eg_data, fit_x, fit_y, slope, 'norg-eg', 'blue')

x_norg_data, norg_t2g_data = np.loadtxt("./norg-t2g-sigma.data", unpack=True, usecols=(0,1))
slope, intercept = fit_near_zero(x_norg_data, norg_t2g_data, threshold)
fit_y = slope * fit_x + intercept
plot_fit(x_norg_data, norg_t2g_data, fit_x, fit_y, slope, 'norg-t2g', 'orange')

x_ctqmc_data, ctqmc_eg_data = np.loadtxt("./ctqmc-eg-sigma.data", unpack=True, usecols=(0,1))
slope, intercept = fit_near_zero(x_ctqmc_data, ctqmc_eg_data, threshold)
fit_y = slope * fit_x + intercept
plot_fit(x_ctqmc_data, ctqmc_eg_data, fit_x, fit_y, slope, 'ctqmc-eg', 'green')

x_ctqmc_data, ctqmc_t2g_data = np.loadtxt("./ctqmc-t2g-sigma.data", unpack=True, usecols=(0,1))
slope, intercept = fit_near_zero(x_ctqmc_data, ctqmc_t2g_data, threshold)
fit_y = slope * fit_x + intercept
plot_fit(x_ctqmc_data, ctqmc_t2g_data, fit_x, fit_y, slope, 'ctqmc-t2g', 'red')
