# Script to make a comparison plot from the RJ-MCMC output
# of the AH and IL methods.

# The script assumes that the current directory contains the AH-output, and the directory with the IL-output is specified.
# We also assume that the comparison is direct - i.e. that all parameters and data are shared between the two models.

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os
import sys

if not os.path.exists('input_file'):
    print('*'*50)
    print('Cannot find file: input_file')
    print('Check that you are running this Python code in the AH outputs directory \n E.g. cd Outputs \n python ../make_plots.py')
    print('*'*50)
    sys.exit(0)


#
if len(sys.argv) != 2:
    print('*'*50)
    print("Syntax: python make_AH_IL_comparison.py <IL_directory>")
    print("where IL_directory is the relative path to the IL-output")
    print('*'*50)
    sys.exit(0)


# Read some basic data that was saved from the input file to "model_data"
# This can be overwritten by either altering this file, or simply hardwiring the various parameters: e.g.
# age_min, age_max = 0, 100

for line in open('input_file','r'):
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

# read in the various data files that were output by the RJ-MCMC script
x, x_err, y, y_err, strat = np.loadtxt('data.dat', unpack=True)
strat = [int(a) for a in strat]
lx, ly = np.loadtxt('credible_lower.dat', unpack=True)
ux, uy = np.loadtxt('credible_upper.dat', unpack=True)
mode_x, mode_y = np.loadtxt('mode.dat', unpack=True)
median_x, median_y = np.loadtxt('median.dat', unpack=True)
av_x, av_y = np.loadtxt('average.dat', unpack=True)
best_x, best_y = np.loadtxt('best_fit.dat', unpack=True)

if not os.path.exists(sys.argv[1]):
    print('Path ' + sys.argv[1] + ' does not exist')
    sys.exit(0)

lx_IL, ly_IL = np.loadtxt(os.path.join(sys.argv[1],'credible_lower.dat'), unpack=True)
ux_IL, uy_IL = np.loadtxt(os.path.join(sys.argv[1],'credible_upper.dat'), unpack=True)
mode_x_IL, mode_y_IL = np.loadtxt(os.path.join(sys.argv[1],'mode.dat'), unpack=True)
median_x_IL, median_y_IL = np.loadtxt(os.path.join(sys.argv[1],'median.dat'), unpack=True)
av_x_IL, av_y_IL = np.loadtxt(os.path.join(sys.argv[1],'average.dat'), unpack=True)
best_x_IL, best_y_IL = np.loadtxt(os.path.join(sys.argv[1],'best_fit.dat'), unpack=True)


print('Building comparative figure...')

fig1, ax1 = plt.subplots(figsize=(14,5))

ax1.fill_between(lx, ly, uy, facecolor='orange', alpha=0.5, label = 'AH %i%% Credible interval' % credible)

ax1.plot(av_x_IL, av_y_IL, 'r', label = 'Average: IL', linewidth=2)
ax1.plot(av_x, av_y, 'b', label = 'Average: AH', linewidth=2)

if 'x_cts_true' in locals():  #see if "true" data are available to plot --- only for synthetic cases.
    plt.plot(x_cts_true,y_cts_true,'k', linewidth=2, label='Real')

ax1.plot(lx_IL, ly_IL, 'g-', label = 'IL %i%% Credible interval' % credible,linewidth=2)
ax1.plot(ux_IL, uy_IL, 'g-',linewidth=2)

ax1.set_xlabel('Age/yr',fontsize=16)
ax1.set_ylabel('Intensity/$\mu$T',fontsize=16)
ax1.xaxis.set_tick_params(labelsize=16)
ax1.yaxis.set_tick_params(labelsize=16)

ax1.set_xlim(age_min, age_max)
ax1.set_ylim(I_min, I_max)
ax1.legend(loc = 'upper right',fontsize=12,labelspacing=0.2)
plt.savefig('AH_IL_comparison.pdf', bbox_inches='tight',pad_inches=0.0)
plt.close(fig1)

