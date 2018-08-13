import numpy as np
from scipy import signal
from scipy.stats import norm
import matplotlib.pyplot as plt
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

# If the RJ-MCMC code has saved some models, load the data

if os.path.exists(output_model_filename):
    output_model = open(output_model_filename,'r')
    discretise_size = int(output_model.readline().split()[0])
    times = [float(x) for x in output_model.readline().split()]
    models = [float(x.split()[0]) for x in output_model.readlines()]
    num_models =  int(len(models) / discretise_size)
    models = np.reshape(models,[num_models,discretise_size])
    output_model.close()
else:
    print('Data for individual models not found - no plot made.')
    sys.exit(0)


# Produces power-spectral-density plots and along with a pdf.
# Processing steps:
# 1. Restricting the range of ages of the models defined over a much greater age range.
# 2. Detrending the resulting time-series
# 3. Applying a butterworth filter (twice) to remove frequencies below highcut and above lowcut.
# 4. Calculate the PSD using the method of Welch, using npad frequencies.

min_age = 500
max_age = 2000

P_edge_max = 50000  #max vertical for plots
P_edge_min = 1000

def process(data):
    fs = 1.0 /( data[1,0] - data[0,0])
    datafilt1=signal.detrend(data, axis=0, type='linear')
    nlint1=datafilt1[:,1]
    npad=2048

    f_nyq = fs / 2. 
    f_lowcut = 1.0 / 400.0
    w_lowcut = f_lowcut / f_nyq
# Begin Phil
    f_highcut = 1.0 / 60.
# End Phil
    w_highcut = f_highcut / f_nyq
# define filter properties
    b,a = signal.butter(2, [w_lowcut,w_highcut], 'band')
# apply filter to signal (twice, to avoid phasing) 
    nlint2 = signal.filtfilt(b,a,nlint1)

    f,P=signal.welch(nlint2,fs,detrend='constant',nperseg=56,nfft=npad)
    T=1/f
    return (T,P)

total_P = []
total_T = []
av_P = np.array(0)

fig3, ax3 = plt.subplots( figsize = (8,5))

plotting_range_min,plotting_range_max = 150, 400
T_grid = np.linspace(plotting_range_min,plotting_range_max,100)

peak_psd_period = np.zeros(num_models)
for i in range(num_models):
    if np.mod(i+1,100) == 0: print('Processing ' + str(i+1) + ' of ' + str(num_models))
    intensity = models[i,:]
    
    restrict_ages = [j for j in range(0,len(times)) if (times[j] > min_age) and (times[j] < max_age) ]
    data = np.zeros([len(restrict_ages),2])
    data[:,0] = np.array(times)[restrict_ages]
    data[:,1] = np.array(intensity)[restrict_ages]
    
    (T,P) = process(data)
# interpolate onto an equally spaced grid for plotting
    P_interp = np.interp(T_grid, T[::-1],P[::-1])
    #ax2.semilogx(T,P,linewidth=0.1,color='black')
    ax3.plot(T,P,linewidth=0.1,color='black',alpha=0.4)
    peak_psd_period[i] = T_grid[np.argmax(P_interp)]
# make a big list of data points which hist2d will compile into a density plot
    total_T.append(T_grid)
    total_P.append(P_interp)
    av_P = av_P + P
    
av_P = av_P / float(num_models)
total_T = np.array(total_T)
total_P = np.array(total_P)

# define the bin edges for the power spectral density:
P_edges = np.linspace(P_edge_min,P_edge_max,100)

fig1, ax1 = plt.subplots(figsize=(8,5))

H = ax1.hist2d(total_T.flatten(),total_P.flatten(), bins=(T_grid, P_edges), normed=True)

cb = plt.colorbar(H[3],ax=ax1)
cb.set_label('Discrete probability',fontsize=16)

ax1.set_ylim([P_edge_min,P_edge_max])
ax1.set_xlabel('Period/yr',fontsize=16)
ax1.set_ylabel('PSD/$\mu$T$^2$y',fontsize=16)
ax1.xaxis.set_tick_params(labelsize=16)
ax1.yaxis.set_tick_params(labelsize=16)
fig1.savefig('PSD_density.pdf', bbox_inches='tight',pad_inches=0.0)
plt.close(fig1)


# add median data:
data=np.loadtxt('median.dat')
data_median = np.array(data)[restrict_ages]
(T_median,P_median) = process(data_median)
spectrum_median_plot, = ax3.plot(T_median,P_median,color='green')
         
# add mode data:
data=np.loadtxt('mode.dat')
data_mode = np.array(data)[restrict_ages]
(T_mode,P_mode) = process(data_mode)
spectrum_mode_plot, = ax3.plot(T_mode,P_mode,color='red')

# add the average model:
data=np.loadtxt('average.dat')
data_av = np.array(data)[restrict_ages]
(T_av,P_av) = process(data_av)
spectrum_average_plot, = ax3.plot(T_av,P_av,color='blue')

print('Peak of average spectra occurs at period %f (yr)'%T_av[np.argmax(P_av)])

# Plot the "average power spectra" over all models read in:
average_spectra_plot, = ax3.plot(T,av_P,linewidth=2,color='orange')

# finish up
ax3.set_ylim([0,P_edge_max])
ax3.set_xlabel('Period/yr',fontsize=16)
ax3.set_ylabel('PSD/$\mu$T$^2$yr',fontsize=16)
ax3.set_xlim([plotting_range_min,plotting_range_max])
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.xaxis.set_tick_params(labelsize=16)
ax3.yaxis.set_tick_params(labelsize=16)
fig3.legend([spectrum_median_plot,spectrum_mode_plot,spectrum_average_plot,average_spectra_plot],['Median','Mode','Average','Average spec'],loc = 'upper right',fontsize=16,labelspacing=0.2)

fig3.savefig('PSD_linear.pdf', bbox_inches='tight',pad_inches=0.0)
plt.close(fig3)

# Make plot of PSD maxima:
fig2, ax2 = plt.subplots( figsize=(8,5))
ax2.hist(peak_psd_period,bins=20,color='blue', edgecolor='white',range=[100,400], density=True,alpha=0.6)

# Fit a normal distribution to the data:
mu, std = norm.fit(peak_psd_period)
print('Best fit normal distribution for maxima for PSD has (mean, std) of ' + str(mu) + ' ' + str(std) )
xmin, xmax = ax2.get_xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
#ax2.plot(x, p, 'k', linewidth=2)
ax2.xaxis.set_tick_params(labelsize=16)
ax2.yaxis.set_tick_params(labelsize=16)
ax2.set_xlabel('Period of PSD maximum/yr',fontsize=16)
ax2.set_ylabel('Discrete probability',fontsize=16)
fig2.savefig('PSD_max.pdf', bbox_inches='tight',pad_inches=0.0)

