# Script to make plots from the RJ-MCMC output

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import numpy as np
from scipy import stats
from matplotlib.colors import LogNorm
import os
import sys

if not os.path.exists('Outputs_Paris700_no_age_errors/input_file'):
    print('*'*50)
    print('Cannot find file: Outputs_Paris700_no_age_errors/input_file')
    print('Check that you are running this Python code in the main directory \n')
    print('*'*50)
    sys.exit(0)


# Read some basic data from the input file
# This can be overwritten by either altering this file, or simply hardwiring the various parameters: e.g.
# age_min, age_max = 0, 100

for line in open('Outputs_Paris700_no_age_errors/input_file','r'):
    if not (line[0] == '#' or line == '\n'): #skip comments or blank lines...
        if line.split()[0].upper() == 'Intensity_prior'.upper():
            I_min, I_max =  float(line.split()[1]),float(line.split()[2])
        if line.split()[0].upper() == 'Age_bounds'.upper():
            age_min, age_max =  float(line.split()[1]),float(line.split()[2])
        if line.split()[0].upper() == 'Num_change_points'.upper():
            K_min, K_max =  int(line.split()[1]), int(line.split()[2])
        if line.split()[0].upper() == 'Credible'.upper():
            credible = float(line.split()[1])
        if line.split()[0].upper() == 'output_model'.upper():
            output_model_filename = line.split()[1]
        if line.split()[0].upper() == 'True_data'.upper():
            true_behaviour_file = line.split()[2]
            x_cts_true,y_cts_true=np.loadtxt(os.path.join(os.pardir,true_behaviour_file),unpack=True)
        if line.split()[0].upper() == 'Plotting_intensity_range'.upper():
            I_min,I_max =  float(line.split()[1]),float(line.split()[2])


plt.figure(figsize=(14,12))
#Part 1:
ax1 = plt.subplot(411)
ax = ax1
dir = 'Outputs_Paris700/'
x, x_err, y, y_err, strat = np.loadtxt(dir+'data.dat', unpack=True)
strat = [int(a) for a in strat]
lx, ly = np.loadtxt(dir+'credible_lower.dat', unpack=True)
ux, uy = np.loadtxt(dir+'credible_upper.dat', unpack=True)
mode_x, mode_y = np.loadtxt(dir+'mode.dat', unpack=True)
median_x, median_y = np.loadtxt(dir+'median.dat', unpack=True)
av_x, av_y = np.loadtxt(dir+'average.dat', unpack=True)

ax.fill_between(lx, ly, uy, facecolor='orange', alpha=0.5, edgecolor='g', label='%i%% credible interval' % credible)
(line, caps, bars) = ax.errorbar(x, y,xerr=x_err, yerr=y_err,fmt='o',color='blue',ecolor='k', elinewidth=1, capthick=0.7, capsize=4, markersize=5)
ax.plot(av_x, av_y, 'r', label = 'Average', linewidth=2)
ax.plot(median_x, median_y, 'purple', linewidth=2, label = 'Median')
ax.plot(mode_x, mode_y, 'blue', linewidth=2, label = 'Mode')

ax.set_ylim(I_min,I_max)
ax.set_xlim(age_min, age_max)
ax.set_title('Posterior distribution: data ages free parameters',fontsize=16)
#ax.set_xlabel('Age/yr',fontsize=16)
ax.set_ylabel('Intensity/$\mu$T',fontsize=16)
plt.setp(line,label="Data")
ax.legend(loc = 'upper right',fontsize=12,labelspacing=0.2)
#ax.xaxis.set_tick_params(labelsize=16)
ax.yaxis.set_tick_params(labelsize=16)
ax.xaxis.grid(True)
plt.setp(ax.get_xticklabels(), visible=False)

#
ax2 = plt.subplot(412,sharex=ax1)
dir = 'Outputs_Paris700_no_age_errors/'
x, x_err, y, y_err, strat = np.loadtxt(dir+'data.dat', unpack=True)
strat = [int(a) for a in strat]
lx, ly = np.loadtxt(dir+'credible_lower.dat', unpack=True)
ux, uy = np.loadtxt(dir+'credible_upper.dat', unpack=True)
mode_x, mode_y = np.loadtxt(dir+'mode.dat', unpack=True)
median_x, median_y = np.loadtxt(dir+'median.dat', unpack=True)
av_x, av_y = np.loadtxt(dir+'average.dat', unpack=True)

ax = ax2
ax.fill_between(lx, ly, uy, facecolor='orange', alpha=0.5, edgecolor='g', label='%i%% credible interval' % credible)
(line, caps, bars) = ax.errorbar(x, y,xerr=x_err, yerr=y_err,fmt='o',color='blue',ecolor='k', elinewidth=1, capthick=0.7, capsize=4, markersize=5)
ax.plot(av_x, av_y, 'r', label = 'Average', linewidth=2)
ax.plot(median_x, median_y, 'purple', linewidth=2, label = 'Median')
ax.plot(mode_x, mode_y, 'blue', linewidth=2, label = 'Mode')

ax.set_ylim(I_min,I_max)
ax.set_xlim(age_min, age_max)
ax.set_title('Posterior distribution: age errors assumed zero',fontsize=16)
#ax.set_xlabel('Age/yr',fontsize=16)
ax.set_ylabel('Intensity/$\mu$T',fontsize=16)
#ax.legend(loc = 'upper right',fontsize=12,labelspacing=0.2)
#ax.xaxis.set_tick_params(labelsize=16)
ax.yaxis.set_tick_params(labelsize=16)
ax.xaxis.grid(True)
plt.setp(ax.get_xticklabels(), visible=False)

#
ax3 = plt.subplot(413,sharex=ax1)
ax = ax3
dir = 'Outputs_Paris700_no_age_errors_twice_F_error/'
x, x_err, y, y_err, strat = np.loadtxt(dir+'data.dat', unpack=True)
strat = [int(a) for a in strat]
lx, ly = np.loadtxt(dir+'credible_lower.dat', unpack=True)
ux, uy = np.loadtxt(dir+'credible_upper.dat', unpack=True)
mode_x, mode_y = np.loadtxt(dir+'mode.dat', unpack=True)
median_x, median_y = np.loadtxt(dir+'median.dat', unpack=True)
av_x, av_y = np.loadtxt(dir+'average.dat', unpack=True)

ax.fill_between(lx, ly, uy, facecolor='orange', alpha=0.5, edgecolor='g', label='%i%% credible interval' % credible)
(line, caps, bars) = ax.errorbar(x, y,xerr=x_err, yerr=y_err,fmt='o',color='blue',ecolor='k', elinewidth=1, capthick=0.7, capsize=4, markersize=5)
ax.plot(av_x, av_y, 'r', label = 'Average', linewidth=2)
ax.plot(median_x, median_y, 'purple', linewidth=2, label = 'Median')
ax.plot(mode_x, mode_y, 'blue', linewidth=2, label = 'Mode')

ax.set_ylim(I_min,I_max)
ax.set_xlim(age_min, age_max)
ax.set_title('Posterior distribution: age errors assumed zero with twice the intensity error budget',fontsize=16)
#ax.set_xlabel('Age/yr',fontsize=16)
ax.set_ylabel('Intensity/$\mu$T',fontsize=16)
#ax.legend(loc = 'upper right',fontsize=12,labelspacing=0.2)
#ax.xaxis.set_tick_params(labelsize=16)
ax.yaxis.set_tick_params(labelsize=16)
ax.xaxis.grid(True)
plt.setp(ax.get_xticklabels(), visible=False)


ax4 = plt.subplot(414,sharex=ax1)
ax = ax4

dir = 'Outputs_Paris700_no_age_errors_min_5muT/'
x, x_err, y, y_err, strat = np.loadtxt(dir+'data.dat', unpack=True)
strat = [int(a) for a in strat]
lx, ly = np.loadtxt(dir+'credible_lower.dat', unpack=True)
ux, uy = np.loadtxt(dir+'credible_upper.dat', unpack=True)
mode_x, mode_y = np.loadtxt(dir+'mode.dat', unpack=True)
median_x, median_y = np.loadtxt(dir+'median.dat', unpack=True)
av_x, av_y = np.loadtxt(dir+'average.dat', unpack=True)

ax.fill_between(lx, ly, uy, facecolor='orange', alpha=0.5, edgecolor='g', label='%i%% credible interval' % credible)
(line, caps, bars) = ax.errorbar(x, y,xerr=x_err, yerr=y_err,fmt='o',color='blue',ecolor='k', elinewidth=1, capthick=0.7, capsize=4, markersize=5)
ax.plot(av_x, av_y, 'r', label = 'Average', linewidth=2)
ax.plot(median_x, median_y, 'purple', linewidth=2, label = 'Median')
ax.plot(mode_x, mode_y, 'blue', linewidth=2, label = 'Mode')

ax.set_ylim(I_min,I_max)
ax.set_xlim(age_min, age_max)
ax.set_title('Posterior distribution: age errors assumed zero and minimum 5 $\mu$T intensity error',fontsize=16)
ax.set_xlabel('Time/yr',fontsize=16)
ax.set_ylabel('Intensity/$\mu$T',fontsize=16)
#ax.legend(loc = 'upper right',fontsize=12,labelspacing=0.2)
#ax.xaxis.set_tick_params(labelsize=16)
#ax.yaxis.set_tick_params(labelsize=16)
ax.xaxis.grid(True)
#plt.setp(ax.get_xticklabels(), visible=False)
# read in the various data files that were output by the RJ-MCMC script


ax.xaxis.set_tick_params(labelsize=16)
ax.yaxis.set_tick_params(labelsize=16)
plt.savefig('Paris_age_sampling_comparison.pdf', bbox_inches='tight',pad_inches=0.4)
plt.close()


