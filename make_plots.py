# Script to make plots from the RJ-MCMC output

import numpy as np
import matplotlib.pyplot as plt
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
        if line.split()[0].upper() == 'Burn_in'.upper():
            Burn_in = int(line.split()[1])
        if line.split()[0].upper() == 'Outputs_directory'.upper():
            outputs_directory = line.split()[1]
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

print('Building plot of data...')
# plot the data with the density binned by using the number of bins of "num_bins"
num_bins = 20

fig1, ax1 = plt.subplots(figsize=(14,6))

unstratified_index = [index for index,item in enumerate(strat) if item == 0]
stratified_index = [index for index,item in enumerate(strat) if item > 0]

if len(unstratified_index) > 0:
    (line, caps, bars) = ax1.errorbar(x[unstratified_index], y[unstratified_index],xerr=x_err[unstratified_index], yerr=y_err[unstratified_index],
     fmt='o',markerfacecolor='blue',markeredgecolor='k', markeredgewidth = 0.6, ecolor='k', elinewidth=1, capsize=4, markersize=7)
    plt.setp(line,label="Unstratified data") #give label to returned line


if len(stratified_index) > 0:
    (line2, caps, bars) = ax1.errorbar(x[stratified_index], y[stratified_index],xerr=x_err[stratified_index], yerr=y_err[stratified_index],
    fmt='o',markerfacecolor='red',markeredgecolor='k', markeredgewidth = 0.6, ecolor='k', elinewidth=1, capsize=4, markersize=7)
    plt.setp(line2,label="Stratified data") #give label to returned line


ax1.set_xlabel('Time/yr',fontsize=16)
ax1.set_ylabel('Intensity/$\mu$T',fontsize=16)
ax1.xaxis.set_tick_params(labelsize=16)
ax1.yaxis.set_tick_params(labelsize=16)

count_colour = 'g'
ax2 = ax1.twinx()
ax2.hist(x,num_bins,alpha=0.2,color=count_colour,edgecolor='white')
ax2.set_ylabel('Count',fontsize=16,color=count_colour)
for tl in ax2.get_yticklabels():
    tl.set_color(count_colour)

ax1.set_xlim(age_min, age_max)
ax2.yaxis.set_tick_params(labelsize=16)

plt.savefig('Data.pdf', bbox_inches='tight',pad_inches=0.0)
plt.close(fig1)



# Make a single plot of the data with mean/mode/median/credible bounds for the posterior
print('Building plot of posterior...')
fig2, ax = plt.subplots (figsize=(14,5))
ax.fill_between(lx, ly, uy, facecolor='orange', alpha=0.5, edgecolor='g', label='%i%% credible interval' % credible)

#a.errorbar(dx[black_pts_index], dy[black_pts_index],xerr=dx_err[black_pts_index], yerr=dn[black_pts_index],fmt='k.', label='Data', elinewidth=0.5)

(line, caps, bars) = ax.errorbar(x, y,xerr=x_err, yerr=y_err,fmt='o',color='blue',ecolor='k', elinewidth=1, capthick=0.7, capsize=4, markersize=5)
plt.setp(line,label="Data") #give label to returned line
    
ax.plot(av_x, av_y, 'r', label = 'Average', linewidth=2)
#ax.plot(best_x, best_y, 'b', linewidth=2, label = 'Best fit')
ax.plot(median_x, median_y, 'purple', linewidth=2, label = 'Median')
ax.plot(mode_x, mode_y, 'blue', linewidth=2, label = 'Mode')
if 'x_cts_true' in locals():  #see if "true" data are available to plot --- only for synthetic cases.
    plt.plot(x_cts_true,y_cts_true,'k', linewidth=2, label='Real')



ax.set_ylim(I_min,I_max)
ax.set_xlim(age_min, age_max)
ax.set_title('Posterior distribution of intensity',fontsize=20)
ax.set_xlabel('Age/yr',fontsize=16)
ax.set_ylabel('Intensity/$\mu$T',fontsize=16)
ax.legend(loc = 'lower right',fontsize=12,labelspacing=0.2)
ax.xaxis.set_tick_params(labelsize=16)
ax.yaxis.set_tick_params(labelsize=16)
plt.savefig('Posterior.pdf', bbox_inches='tight',pad_inches=0.4)
plt.close(fig2)


# Make a plot of the histogram of the number of change points
print('Building plot of change points...')
fig3, ax = plt.subplots (figsize=(8,5))
k_count = k_count/np.sum(k_count) #normalise
ax.bar(k_index,k_count,align='center')
#ax.set_xticks(k_index[::2])

ax.set_title('Vertices Histogram',fontsize=16)
ax.set_xlabel('Number of vertices',fontsize=16)
ax.set_ylabel('Discrete probability',fontsize=16)
ax.xaxis.set_tick_params(labelsize=16)
ax.yaxis.set_tick_params(labelsize=16)
plt.savefig('K_histogram.pdf', bbox_inches='tight',pad_inches=0.4)
plt.close(fig3)


# Make a plot of the age of the change points
num_bins = 500
vertices = np.loadtxt('changepoints.dat')
fig4, ax = plt.subplots (figsize=(14,3))
ax.hist(vertices, bins = num_bins)
ax.set_title('Vertex position Histogram',fontsize=20)
ax.set_xlabel('Age/yr',fontsize=16)
ax.set_ylabel('Count',fontsize=16)
ax.set_xlim(age_min, age_max)
ax.xaxis.set_tick_params(labelsize=16)
ax.yaxis.set_tick_params(labelsize=16)
plt.savefig('Change_point_histogram.pdf', bbox_inches='tight',pad_inches=0.4)
plt.close(fig4)

# Make a plot of the misfit
print('Building plot of misfit...')

iterations, misfit = np.loadtxt('misfit.dat',unpack=True)
fig5, ax = plt.subplots (figsize=(8,5) )
ax.plot(iterations, misfit,'k')
ax.set_yscale('log')
ax.set_title('Misfit against iteration count',fontsize=16)
ax.set_xlabel('Iteration count',fontsize=16)
ax.set_ylabel('Misfit',fontsize=16)
ax.xaxis.set_tick_params(labelsize=16)
ax.yaxis.set_tick_params(labelsize=16)

# add red bar to indicate the burn-in end
ax.bar(Burn_in,height=misfit.max(),width=iterations.max()/100,bottom = 0, align = 'center',color='red')
plt.savefig('Misfit.pdf', bbox_inches='tight',pad_inches=0.4)
plt.close(fig5)

# Make a plot of the density
threshold = 0.01
print('Building plot of density...')
fig6, ax = plt.subplots ( figsize=(14,5))

ax.set_title('Intensity density',fontsize=18)
ax.set_ylabel('Intensity/$\mu$T')
#f = open('intensity_density.dat', 'r')
#discretise_size, NBINS = [int(x) for x in f.readline().split()]
#density_data =  [list(map(float, x.split())) for x in f.readlines()]
#f.close()
#x_density,y_density,intensity_density = list(zip(*density_data))
#int_density = np.reshape(intensity_density,[discretise_size,NBINS])
#x_density = np.reshape(x_density,[discretise_size,NBINS])
#y_density = np.reshape(y_density,[discretise_size,NBINS])
#int_density = np.transpose(int_density)
#plt.imshow(int_density, origin='lower',cmap = cm.jet,extent=(x_density[0,0],x_density[-1,0],y_density[0,0],y_density[0,-1]), aspect='auto')

f = open('intensity_density.dat', 'r')
discretise_size, NBINS = [int(x) for x in f.readline().split()]
density_data =  [list(map(float, x.split())) for x in f.readlines()]
f.close()
x_density,y_density,intensity_density = list(zip(*density_data))
int_density = np.reshape(intensity_density,[discretise_size,NBINS])
x_density = np.reshape(x_density,[discretise_size,NBINS])
y_density = np.reshape(y_density,[discretise_size,NBINS])

int_density = np.transpose(int_density)
x_density = np.transpose(x_density)
y_density = np.transpose(y_density)
int_density_refined = int_density.copy()

int_density_refined[ int_density_refined < threshold] = 0.0

plt.imshow(int_density_refined, origin='lower',cmap = cm.jet,extent=(x_density[0,0],x_density[0,-1],y_density[0,0],y_density[-1,0]), aspect='auto', norm=LogNorm(vmin=0.01, vmax=0.6),interpolation="nearest")
#ax2.set_ylabel('Intensity/$\mu$T',fontsize=16)



#plt.xlim(x_density[0,0],x_density[-1,0])
#plt.ylim(y_density[0,0],y_density[0,-1])
plt.xlabel('Age/yr',fontsize=16)
plt.ylabel('Intensity/$\mu$T',fontsize=16)

cb = plt.colorbar(ticks=[0.01, 0.1, 0.2, 0.3, 0.4, 0.5, .6],
                  orientation='vertical', format='$%.2f$')

    #cb = plt.colorbar(ticks=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6,0.7, 0.8,0.9, 1],
#             orientation='vertical')
cb.set_label('Probability',fontsize=16)
#plt.clim(0, 1)
ax.xaxis.set_tick_params(labelsize=16)
ax.yaxis.set_tick_params(labelsize=16)
#plt.plot(x_orig,y_orig,'white', linewidth=1, label='Real')
plt.savefig('density.pdf', bbox_inches='tight',pad_inches=0.4)
plt.close(fig6)

# Make a plot of the intensity error - to check against the prior assumption
print('Building plot of intensity error...')

fig, ax = plt.subplots ( figsize=(8,5) )
# interpolate the average curve
interp_data = np.interp(x,av_x,av_y)
weighted_errors = (interp_data - y) / y_err

n, bins, patches = plt.hist(weighted_errors,bins=100,density=1,edgecolor='w',range=[-5,5])

x_smooth = np.linspace(-5,5,1000)
normal_distribution = stats.norm.pdf(x_smooth, 0, 1)
ax.plot(x_smooth,normal_distribution,'r')

ax.set_title('Weighted intensity error',fontsize=16)
ax.set_xlabel('$\sigma^{-1} \; \Delta F$',fontsize=16)
ax.set_ylabel('Discrete probability', fontsize=16)
ax.xaxis.set_tick_params(labelsize=16)
ax.yaxis.set_tick_params(labelsize=16)
ax.set_xlim(-5,5)

plt.savefig('histogram_weighted_errors.pdf', bbox_inches='tight',pad_inches=0.4)
plt.close(fig)




# Make a 3-part joint plot
intensity_range = 70
threshold = 0.01
print('Building composite plot...')
plt.figure(figsize=(14,15))

#Part 1:
ax1 = plt.subplot(311)
ax1.fill_between(lx, ly, uy, facecolor='orange', alpha=0.5, edgecolor='g', label='%i%% credible interval' % credible)

unstratified_index = [index for index,item in enumerate(strat) if item == 0]
stratified_index = [index for index,item in enumerate(strat) if item > 0]

if len(unstratified_index) > 0:
    (line, caps, bars) = ax1.errorbar(x[unstratified_index], y[unstratified_index],xerr=x_err[unstratified_index], yerr=y_err[unstratified_index],
    fmt='o',color='blue',ecolor='k', elinewidth=1, capthick=0.7, capsize=4, markersize=5)
    plt.setp(line,label="Unstratified data") #give label to returned line

if len(stratified_index) > 0:
    (line2, caps, bars) = ax1.errorbar(x[stratified_index], y[stratified_index],xerr=x_err[stratified_index], yerr=y_err[stratified_index],
    fmt='o',color='red',ecolor='k', elinewidth=1, capthick=0.7, capsize=4, markersize=5)
    plt.setp(line2,label="Stratified data") #give label to returned line
#for cap in caps:
#    cap.set_markeredgewidth(0.5)

ax1.plot(av_x, av_y, 'r', label = 'Average', linewidth=2)
ax1.plot(median_x, median_y, 'purple', linewidth=2, label = 'Median')
ax1.plot(mode_x, mode_y, 'blue', linewidth=2, label = 'Mode')
if 'x_cts_true' in locals():  #see if "true" data are available to plot --- only for synthetic cases.
    ax1.plot(x_cts_true,y_cts_true,'k', linewidth=2, label='Real')

ax1.set_ylim(I_min,I_max)
ax1.set_title('Posterior distribution of intensity',fontsize=20)
#ax1.set_xlabel('Time/yr',fontsize=16)
ax1.set_ylabel('Intensity/$\mu$T',fontsize=16)
ax1.legend(numpoints =1, loc = 'upper right',fontsize=12,labelspacing=0.2)
ax1.yaxis.set_tick_params(labelsize=16)
# Make x-tick labels invisible.
plt.setp(ax1.get_xticklabels(), visible=False)


# Part 2
ax2 = plt.subplot(312,sharex=ax1)
ax2.set_title('Intensity density',fontsize=20)

f = open('intensity_density.dat', 'r')
discretise_size, NBINS = [int(x) for x in f.readline().split()]
density_data =  [list(map(float, x.split())) for x in f.readlines()]
f.close()
x_density,y_density,intensity_density = list(zip(*density_data))
int_density = np.reshape(intensity_density,[discretise_size,NBINS])
x_density = np.reshape(x_density,[discretise_size,NBINS])
y_density = np.reshape(y_density,[discretise_size,NBINS])

int_density = np.transpose(int_density)
x_density = np.transpose(x_density)
y_density = np.transpose(y_density)
int_density_refined = int_density.copy()

int_density_refined[ int_density_refined < threshold] = 0.0

plt.imshow(int_density_refined, origin='lower',cmap = cm.jet,extent=(x_density[0,0],x_density[0,-1],y_density[0,0],y_density[-1,0]), aspect='auto', norm=LogNorm(vmin=0.01, vmax=0.6),interpolation="nearest")
ax2.set_ylabel('Intensity/$\mu$T',fontsize=16)


ax2.xaxis.set_tick_params(labelsize=16)
ax2.yaxis.set_tick_params(labelsize=16)
if 'x_cts_true' in locals():  #see if "true" data are available to plot --- only for synthetic cases.
    ax2.plot(x_cts_true,y_cts_true,'k', linewidth=2, label='Real')
ax2.set_ylim(I_min,I_max)

# Part 3:
num_bins = 200
ax3 = plt.subplot(313,sharex=ax1)
ax3.hist(vertices, bins = num_bins, histtype='bar', edgecolor='white', linewidth=0.1, color='b')
ax3.set_ylabel('Count',fontsize=16,color='b')
ax3.yaxis.set_tick_params(labelsize=16,colors='b')
# Make x-tick labels invisible.
plt.setp(ax2.get_xticklabels(), visible=False)
part3_title = 'Vertex position'

# KL divergence
# check to see if prior distribution exists
priors_directory = os.path.join(os.pardir,outputs_directory+'_prior_sampling')
if os.path.exists(priors_directory):
    part3_title = part3_title + '/Kullback-Leibler divergence'

# load priors distribution
    f = open(priors_directory+'/intensity_density.dat', 'r')
    discretise_size, NBINS = [int(x) for x in f.readline().split()]
    density_data =  [list(map(float, x.split())) for x in f.readlines()]
    f.close()
    x_density_prior,y_density_prior,intensity_density_prior = list(zip(*density_data))
    int_density_prior = np.reshape(intensity_density_prior,[discretise_size,NBINS])
    x_density_prior = np.reshape(x_density_prior,[discretise_size,NBINS])
    y_density_prior = np.reshape(y_density_prior,[discretise_size,NBINS])
    int_density_prior = np.transpose(int_density_prior)
    x_density_prior = np.transpose(x_density_prior)
    y_density_prior = np.transpose(y_density_prior)
    if not np.array_equal(y_density_prior, y_density) or not np.array_equal(x_density_prior,x_density):
        print('Prior distribution exists but x and y ranges differ.')
        sys.exit(0)
              
# For each time, simply compare the posterior pdf (which integrates to 1) and the prior pdf (which integrates to 1).
    KL = np.zeros( discretise_size )
    for i in range(0, discretise_size):
        int_density_prior[ int_density_prior == 0] = 1e-8  #replace zero values with a small value to avoid dividing by zero.
        KL[i] = stats.entropy(int_density[:,i], qk=int_density_prior[:,i])

    ax3_twin = ax3.twinx()
    ax3_twin.plot(x_density[1,:], KL, 'r')
    ax3_twin.set_ylabel('KL divergence', fontsize=16, color='r')
    ax3_twin.tick_params('y', colors='r',labelsize=16)
    ax3_twin.tick_params(direction='out', pad=5)

# produce an image of the prior intensity:
    threshold_prior = 0.001
    int_density_prior_refined = int_density_prior.copy()
    int_density_prior_refined[ int_density_prior_refined < threshold_prior] = 0.0
    fig8, ax8 = plt.subplots ( figsize=(8,5) )
    plt.sca(ax8)
    plt.imshow(int_density_prior_refined, origin='lower',cmap = cm.jet,extent=(x_density_prior[0,0],x_density_prior[0,-1],y_density_prior[0,0],y_density_prior[-1,0]), aspect='auto', norm=LogNorm(vmin=0.001, vmax=0.1),interpolation="nearest")
    ax8.set_xlim(x_density_prior[0,0],x_density_prior[0,-1])
    ax8.set_ylim(y_density_prior[0,0],y_density_prior[-1,0])
    ax8.set_xlabel('Age/yr',fontsize=16)
    ax8.set_ylabel('Intensity/$\mu$T',fontsize=16)

    cb = plt.colorbar(ticks=[0.001, 0.002, 0.003, 0.004, 0.005, 0.006,0.007, 0.008,0.009, 0.01],orientation='vertical')
    cb.set_label('Probability',fontsize=16)
    ax8.xaxis.set_tick_params(labelsize=16)
    ax8.yaxis.set_tick_params(labelsize=16)
    plt.savefig('Prior_density.pdf', bbox_inches='tight',pad_inches=0.4)
    plt.close(fig8)
    plt.sca(ax1)
    # The x-axis data keeps being updated with new curves - define at the end the x-range:
    
else:
    print('KL divergence ignored as no prior distribution found')
            
ax3.set_xlim(age_min, age_max)
ax3.xaxis.set_tick_params(labelsize=16)
ax3.set_xlabel('Age/yr',fontsize=16)
ax3.xaxis.grid(True)
ax3.tick_params(top = False)
ax1.xaxis.grid(True)
ax2.xaxis.grid(True)
ax1.tick_params(top = False)
ax2.tick_params(top = False)

# Adjust plot sizes to place colour bar
plt.subplots_adjust(bottom=0.05, right=0.9, top=0.95)
cax = plt.axes([0.93, 0.37, 0.03, 0.265])
# The numbers in the square brackets of add_axes refer to [left, bottom, width, height], where the coordinates are just fractions that go from 0 to 1 of the plotting area.
cb = plt.colorbar(ticks=[0.01, 0.1, 0.2, 0.3, 0.4, 0.5, .6],
                  orientation='vertical', format='$%.2f$',cax=cax)
cb.ax.tick_params(labelsize=16)
ax3.set_title(part3_title,fontsize=20)

plt.savefig('joint_plot.pdf', bbox_inches='tight',pad_inches=0.0)
plt.close()

# If the RJ-MCMC code has saved some models, make a plot of these
import os

if os.path.exists(output_model_filename):
    output_model = open(output_model_filename,'r')
    discretise_size = int(output_model.readline().split()[0])
    x_ages = [float(x) for x in output_model.readline().split()]
    models = [float(x.split()[0]) for x in output_model.readlines()]
    num_models =  int(len(models) / discretise_size)
    models = np.reshape(models,[num_models,discretise_size])
    output_model.close()
    print('Building plot of ' + str(num_models) + ' individual models...')
    fig, ax1 = plt.subplots ( figsize=(14,5))
    for i in range(0, num_models):
        ax1.plot( x_ages,models[i,:], 'r',alpha=0.1)
    ax1.plot( av_x, av_y, 'k',  linewidth=2)
    ax1.tick_params(labelsize=16)
    ax1.set_xlabel('Age/yr',fontsize=16)
    ax1.set_ylabel('Intensity/$\mu$T',fontsize=16)
    ax1.set_xlim([age_min,age_max])
    ax1.set_ylim([I_min,I_max])
#ax1.tick_params(top = False)

    plt.savefig('individual_models.pdf', bbox_inches='tight',pad_inches=0.0)
    plt.close(fig)
else:
    print('Data for individual models not found - no plot made.')


# Make dF/dt plots...
print('Building plot of dF/dt...')

lx, ly = np.loadtxt('credible_dFdt_lower.dat', unpack=True)
ux, uy = np.loadtxt('credible_dFdt_upper.dat', unpack=True)
mode_x, mode_y = np.loadtxt('mode_dFdt.dat', unpack=True)
median_x, median_y = np.loadtxt('median_dFdt.dat', unpack=True)
av_x, av_y = np.loadtxt('average_dFdt.dat', unpack=True)


fig1, ax1 = plt.subplots(figsize=(14,6))
ax1.fill_between(lx, ly, uy, facecolor='orange', alpha=0.5, edgecolor='g', label='%i%% credible interval' % credible)

ax1.plot(av_x, av_y, 'r', label = 'Average', linewidth=2)
#ax.plot(best_x, best_y, 'b', linewidth=2, label = 'Best fit')
#ax1.plot(median_x, median_y, 'purple', linewidth=2, label = 'Median')
#ax1.plot(mode_x, mode_y, 'blue', linewidth=2, label = 'Mode')

ax1.set_xlabel('Date/yr',fontsize=16)
ax1.set_ylabel('Rate of change/$\mu T yr^{-1}$',fontsize=16)
ax1.xaxis.set_tick_params(labelsize=16)
ax1.yaxis.set_tick_params(labelsize=16)
ax1.legend(loc = 'lower right',fontsize=12,labelspacing=0.2)
#min_dFdt = [min(np.abs(ly[i]),np.abs(uy[i])) if uy[i] * ly[i] > 0 else np.NaN for i in range(0,len(lx))]
#ax1.plot(lx,min_dFdt,color='red')
ax1.set_ylim([-2,2])
plt.savefig('dFdt.pdf', bbox_inches='tight',pad_inches=0.0)
plt.close(fig1)
