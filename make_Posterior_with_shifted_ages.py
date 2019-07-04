# Script to make a Posterior plot with shifted ages from the RJ-MCMC output
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
from scipy import stats
from matplotlib.colors import LogNorm
import os
import sys

if not os.path.exists('input_file'):
    print('*'*50)
    print('Cannot find file: input_file')
    print('Check that you are running this Python code in the outputs directory \n E.g. cd Outputs \n python ../make_plots.py')
    print('*'*50)
    sys.exit(0)


# Read some basic data from the input file
# This can be overwritten by either altering this file, or simply hardwiring the various parameters: e.g.
# age_min, age_max = 0, 100

for line in open('input_file','r'):
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
# overwrite age limits by plotting information if given in the inputfile
        if line.split()[0].upper() == 'Plotting_age_range'.upper():
            age_min, age_max =  float(line.split()[1]),float(line.split()[2])
        if line.split()[0].upper() == 'Burn_in'.upper():
            Burn_in = int(line.split()[1])
# read in the various data files that were output by the RJ-MCMC script

x, x_err, y, y_err, strat = np.loadtxt('data.dat', unpack=True)
strat = [int(a) for a in strat]
lx, ly = np.loadtxt('credible_lower.dat', unpack=True)
ux, uy = np.loadtxt('credible_upper.dat', unpack=True)
mode_x, mode_y = np.loadtxt('mode.dat', unpack=True)
median_x, median_y = np.loadtxt('median.dat', unpack=True)
av_x, av_y = np.loadtxt('average.dat', unpack=True)
best_x, best_y = np.loadtxt('best_fit.dat', unpack=True)
k_index, k_count = np.loadtxt('k_histogram.dat',unpack=True)


# Make a single plot of the data with mean/mode/median/credible bounds for the posterior
print('Building plot of posterior...')
fig2, ax = plt.subplots (figsize=(14,5))
ax.fill_between(lx, ly, uy, facecolor='orange', alpha=0.5, edgecolor='g', label='%i%% credible interval' % credible)

# Calculate average of sample age/intensity
new_age, new_intensity = np.zeros(len(x)), np.zeros(len(x))
for index in range(len(x)):
    dist = np.loadtxt('Joint_distribution_data/Sample_'+str('{:04}'.format(1+index))+'.dat')
    new_age[index], new_intensity[index] = np.mean(dist[:,0]), np.mean(dist[:,1])

(line, caps, bars) = ax.errorbar(x, y,xerr=x_err, yerr=y_err,fmt='o',color='black',ecolor='k', elinewidth=1, capthick=0.7, capsize=4, markersize=5, alpha=0.3)

plt.setp(line,label="Original data") #give label to returned line

ax.plot( new_age, new_intensity, marker = 'o', linestyle = '',color='blue',label='Shifted data')

ax.plot(av_x, av_y, 'r', label = 'Average', linewidth=2)
#ax.plot(best_x, best_y, 'b', linewidth=2, label = 'Best fit')
#ax.plot(median_x, median_y, 'purple', linewidth=2, label = 'Median')
#ax.plot(mode_x, mode_y, 'blue', linewidth=2, label = 'Mode')
if 'x_cts_true' in locals():  #see if "true" data are available to plot --- only for synthetic cases.
    plt.plot(x_cts_true,y_cts_true,'k', linewidth=2, label='Real')

ax.set_ylim(I_min,I_max)
ax.set_xlim(age_min, age_max)
ax.set_title('Posterior intensity with shifted data',fontsize=20)
ax.set_xlabel('Time/yr',fontsize=16)
ax.set_ylabel('Intensity/$\mu$T',fontsize=16)
ax.legend(loc = 'upper left',fontsize=12,labelspacing=0.2)
ax.xaxis.set_tick_params(labelsize=16)
ax.yaxis.set_tick_params(labelsize=16)
plt.savefig('Posterior_with_shifted_data.pdf', bbox_inches='tight',pad_inches=0.4)
plt.close(fig2)

f= open("Posterior_ages.txt","w+")
for index in range(len(x)):
	f.write("%f %f\n" %(new_age[index],new_intensity[index]))
f.close()


