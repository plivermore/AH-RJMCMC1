import numpy as np
import pandas as pd
import seaborn as sns; sns.set(style="white", color_codes=True)
import matplotlib.pyplot as plt
import os
import sys


if not os.path.exists('input_file'):
    print('*'*30)
    print('Run this Python code in the outputs directory \n E.g. python ../make_joint_distribution_plots.py')
    print('*'*30)
    sys.exit(0)

batch_col = 0

for line in open('input_file'):
    if not (line[0] == '#' or line == '\n'):
        if line.split()[0].upper() == 'Data_file'.upper():
            Data_filename = line.split()[1].rstrip();
        if line.split()[0].upper() == 'File_format'.upper():
            id, age_col, d_age_col, F_col, dF_col, Strat_col = int(line.split()[1]),int(line.split()[2]),int(line.split()[3]),int(line.split()[4]), int(line.split()[5]), int(line.split()[6])
        if line.split()[0].upper() == 'Intensity_prior'.upper():
            I_min,I_max =  float(line.split()[1]),float(line.split()[2])
        if line.split()[0].upper() == 'True_data'.upper():
            true_behaviour_file = line.split()[1]
        if line.split()[0].upper() == 'Plotting_intensity_range'.upper():
            I_min,I_max =  float(line.split()[1]),float(line.split()[2])
        if line.split()[0].upper() == 'Batch_generate_joint_distributions'.upper():
            batch_col =  int(line.split()[1])

if batch_col == 0:
    print('Batch column not specified in inputfile')
    sys.exit(0)

# Here is the filename of the (noisy) data file
filename = os.path.join(os.pardir,Data_filename)

# Number of bins for plotting
num_bins=30

# load the data
data = np.loadtxt(filename,usecols=(age_col, d_age_col, F_col, dF_col, batch_col), unpack=False, comments='#')
label = np.loadtxt(filename,usecols=(id), unpack=False, comments='#',dtype=np.str)


for index in range(0,data.shape[0]):
    if data[index,4] == 1:  #make a plot
        print("Making joint plot for sample " + str(index))

        noisy_pt = [data[index-1,0], data[index-1,2]]
        errs = [data[index-1,1], data[index-1,3]]

# load the (true) data file if available
        if 'true_behaviour_file' in locals():
            data = np.loadtxt(os.path.join(os.pardir,true_behaviour_file), usecols=(0,2))
            true_pt = [data[index-1,0], data[index-1,1]]

# load the relevant joint distribution file
        dist = np.loadtxt('Joint_distribution_data/Sample_'+str('{:04}'.format(index))+'.dat')
        joint_data = pd.DataFrame(data=dist,columns=["age","intensity"])

# create plot
        sns.set_context("paper", font_scale=1.8)
        sns.set_style({"axes.linewidth": 0.7})
        sns.set_style("ticks")
        g = sns.JointGrid(x='age',y='intensity',data=joint_data)
        g.ax_marg_x.hist(dist[:,0], bins=np.linspace(noisy_pt[0]-errs[0],noisy_pt[0]+errs[0],num_bins),alpha=0.5,color='g')
        (n, bins, patches) = g.ax_marg_y.hist(dist[:,1], alpha=.5, orientation="horizontal", bins=num_bins,color='g',density=True)
        g.plot_joint(plt.hexbin, gridsize=num_bins,  cmap="Greens", mincnt = 1)

        if 'true_behaviour_file' in locals(): g.ax_joint.plot(true_pt[0],true_pt[1],'^',color='red',markersize=10,label='True value')
# plot (noisy) data point with errors
        (line, caps, bars) = g.ax_joint.errorbar(noisy_pt[0],noisy_pt[1],xerr=errs[0],yerr=errs[1],color='blueviolet',fmt='s',label=label[index],capsize=5,elinewidth=2,capthick=0.7)

#plot histogram of age prior:
        g.ax_marg_x.bar(noisy_pt[0]-errs[0],height=np.shape(dist)[0]/num_bins,width=2.*errs[0],alpha=0.5,color='darkorange',align='edge')

        y_plot_range = max(errs[1]+1,7)
        x_plot_range = errs[0]+5

        g.ax_joint.set_ylim([min(dist[:,1].min()-2,noisy_pt[1]-y_plot_range),max(dist[:,1].max()+2,noisy_pt[1]+y_plot_range)])
        g.ax_joint.set_xlim([noisy_pt[0]-x_plot_range,noisy_pt[0]+x_plot_range])

#plot histogram of intensity 'prior':
        sigma = errs[1]
        mu = noisy_pt[1]
        q = np.linspace(mu - 4 * sigma, mu + 4 * sigma,1000)
        curve = 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (q - mu)**2 / (2 * sigma**2))

        g.ax_marg_y.fill_betweenx(q,curve,0,alpha=0.5,color='darkorange')


        print('Min/Max value of age from joint distribution are: ' + str(dist[:,0].min()) + ' ' + str(dist[:,0].max()) )
        print('Min/Max value of age from Prior distribution are: ' + str(noisy_pt[0]-errs[0]) + ' ' + str(noisy_pt[0]+errs[0]) )

        g.ax_joint.set_ylabel('Intensity/$\mu$T')
        g.ax_joint.set_xlabel('Time/yr')
        g.ax_joint.legend(loc='lower center')
        g.savefig('Joint_'+label[index] + '_hex.pdf')

