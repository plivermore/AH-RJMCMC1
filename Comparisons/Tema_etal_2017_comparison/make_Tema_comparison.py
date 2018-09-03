#!/usr/bin/env python
# Script to make a comparison plot for RJMCMC/Tema with the analysis in Tema et al. 2017

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import numpy as np
from scipy import stats
from matplotlib.colors import LogNorm
import os
import sys
RJMCMC_source = '../../Outputs_Tema_etal_2017_dataset/'

# Read some basic data from the input file
# This can be overwritten by either altering this file, or simply hardwiring the various parameters: e.g.
# age_min, age_max = 0, 100

for line in open(RJMCMC_source+'input_file','r'):
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
        if line.split()[0].upper() == 'Plotting_intensity_range'.upper():
            I_min,I_max =  float(line.split()[1]),float(line.split()[2])
        if line.split()[0].upper() == 'Burn_in'.upper():
            Burn_in = int(line.split()[1])
# read in the various data files that were output by the RJ-MCMC script

x, x_err, y, y_err, strat = np.loadtxt(RJMCMC_source+'data.dat', unpack=True)
strat = [int(a) for a in strat]
lx, ly = np.loadtxt(RJMCMC_source+'credible_lower.dat', unpack=True)
ux, uy = np.loadtxt(RJMCMC_source+'credible_upper.dat', unpack=True)
mode_x, mode_y = np.loadtxt(RJMCMC_source+'mode.dat', unpack=True)
median_x, median_y = np.loadtxt(RJMCMC_source+'median.dat', unpack=True)
av_x, av_y = np.loadtxt(RJMCMC_source+'average.dat', unpack=True)
best_x, best_y = np.loadtxt(RJMCMC_source+'best_fit.dat', unpack=True)

# Make a single plot of the data with mean/mode/median/credible bounds for the posterior
print('Building comparison plot ...')
fig2, ax = plt.subplots (figsize=(14,5))

ax.fill_between(lx, ly, uy, facecolor='orange', alpha=0.5, edgecolor='g', label='%i%% AH-RJMCMC credible interval' % credible)

#a.errorbar(dx[black_pts_index], dy[black_pts_index],xerr=dx_err[black_pts_index], yerr=dn[black_pts_index],fmt='k.', label='Data', elinewidth=0.5)

(line, caps, bars) = ax.errorbar(x, y,xerr=x_err, yerr=y_err,fmt='o',color='k',ecolor='k', elinewidth=1, capthick=0.7, capsize=4, markersize=5)
plt.setp(line,label="Data") #give label to returned line
ax.plot(av_x, av_y, 'r', label = 'Average AH-RJMCMC posterior', linewidth=2)
ax.plot(median_x, median_y, 'g', label = 'Median AH-RJMCMC posterior', linewidth=2)
#ax.plot(best_x, best_y, 'b', linewidth=2, label = 'Best fit')
#ax.plot(median_x, median_y, 'purple', linewidth=2, label = 'Median')
#ax.plot(mode_x, mode_y, 'blue', linewidth=2, label = 'Mode')

# Load curve from Tema et al (2017):
x_Tema, y_Tema, ly_Tema, uy_Tema  = np.loadtxt('Tema_2017_curve.txt', unpack=True, usecols=(0,7,8,9) )

ax.fill_between(x_Tema, ly_Tema, uy_Tema, facecolor='grey', alpha=0.5, edgecolor='g', label=r'95% Tema et al. (2017) credible interval')
ax.plot(x_Tema, y_Tema, 'b', label = r'Average Tema et al. (2017) posterior', linewidth=2)

# add our own Hawaiian curve
x_Hawaii, y_Hawaii  = np.loadtxt('../../Outputs_Hawaii/average.dat', unpack=True )
#ax.plot(x_Hawaii, y_Hawaii, 'skyblue', label = r'Average AH-RJMCMC Hawaii', linewidth=2)

ax.set_ylim(I_min,I_max)
ax.set_xlim(age_min, age_max)
ax.set_title(r'Time dependence of intensity: dataset of Tema et al. (2017)',fontsize=20)
ax.set_xlabel('Time/yr',fontsize=16)
ax.set_ylabel('Intensity/$\mu$T',fontsize=16)
ax.legend(loc = 'upper left',fontsize=12,labelspacing=0.2)
ax.xaxis.set_tick_params(labelsize=16)
ax.yaxis.set_tick_params(labelsize=16)
plt.savefig('Tema_dataset_comparison.pdf', bbox_inches='tight',pad_inches=0.4)
plt.close(fig2)


