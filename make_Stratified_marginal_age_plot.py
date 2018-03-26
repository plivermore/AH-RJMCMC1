import numpy as np
import matplotlib.pyplot as plt

offset_plot = 0
fig, ax = plt.subplots( figsize=(6,6))
Num_bins = 20
colors = ['red', 'blue', 'brown', 'plum', 'purple', 'grey', 'darkgreen', 'cyan', 'black', 'skyblue']
i = 0

fig2, ax2 = plt.subplots( figsize=(6,6))

age, intensity = np.loadtxt('data.dat',unpack=True ,usecols=(0,2))
samples_order = [i for i in range(len(age)) if age[i] == 1665 ]


print('Making age histogram...')
print('Number of age samples is ' + str(len(samples_order)) )

# plot samples
for index in range(len(samples_order)):
    #print( 1+samples_order[index] )
    dist = np.loadtxt('Joint_distribution_data/Sample_'+str('{:04}'.format(1+samples_order[index]))+'.dat')
    (pts,edges) = np.histogram(dist[:,0],density=True,bins=Num_bins)
    midpoints = edges[:-1] + np.diff(edges)/2
    #ax.plot([1550,1800],[offset_plot,offset_plot],color=colors[index])
    # plot prior
    
    ax.fill_between([1665-85,1665+85],offset_plot, offset_plot+5 * 1/(2 * 85.),color=colors[index])
    ax.fill_between(midpoints,offset_plot + pts*5, offset_plot,color=colors[index])
    #ax.plot(midpoints,offset_plot + pts*10,color=colors[index-12])
    #ax.plot()
# add to the plot of "moved" ages:
    ax2.plot( np.mean( dist[:,0] ), np.mean( dist[:,1] ), 's', color=colors[index])
    #print( np.mean( dist[:,0] ), np.mean( dist[:,1] ) )
    offset_plot += 0.2

(ymin, ymax) = ax.get_ylim()
ax.set_ylim([ymin,ymax])
ax.plot([1665,1665],[ymin-0.1, ymax+0.1],'k--')
ax.set_ylabel('')
ax.set_yticklabels([])
ax.tick_params(left = 'off')
ax.xaxis.set_tick_params(labelsize=16)
ax.set_ylabel('Stratified sample index',fontsize=16)
ax.set_xlabel('Age/yr',fontsize=16)
x_lims = ax.get_xlim()
fig.savefig('marginal_ages.pdf')
plt.close(fig)


# Make a single plot of the average and old estimates of the ages
x, x_err, y, y_err, strat = np.loadtxt('data.dat', unpack=True)
strat = [int(a) for a in strat]
lx, ly = np.loadtxt('credible_lower.dat', unpack=True)
ux, uy = np.loadtxt('credible_upper.dat', unpack=True)
mode_x, mode_y = np.loadtxt('mode.dat', unpack=True)
median_x, median_y = np.loadtxt('median.dat', unpack=True)
av_x, av_y = np.loadtxt('average.dat', unpack=True)

ax2.fill_between(lx, ly, uy, facecolor='orange', alpha=0.5, edgecolor='g', label='credible interval')

for index in range(len(samples_order)):
    index2 = samples_order[index]
    (line2, caps, bars) = ax2.errorbar(x[index2], y[index2],xerr=x_err[index2], yerr=y_err[index2],
        fmt='o',color=colors[index],ecolor=colors[index], elinewidth=1, capthick=0.7, capsize=4, markersize=5)
#plt.setp(line2,label="Stratified data") #give label to returned line

ax2.plot(av_x, av_y, 'r', label = 'Average Posterior', linewidth=2)
I_min,I_max = 35, 65
ax2.set_ylim([I_min,I_max])
ax2.set_xlim(x_lims)
ax2.set_ylabel('Intensity/$\mu$T',fontsize=16)
ax2.legend(numpoints =1, loc = 'upper right',fontsize=16,labelspacing=0.2)
ax2.yaxis.set_tick_params(labelsize=16)
ax2.xaxis.set_tick_params(labelsize=16)
# Make x-tick labels invisible.
#plt.setp(ax2.get_xticklabels(), visible=False)
ax2.set_xlabel('Age/yr',fontsize=16)
fig2.savefig('age_comparison.pdf')
plt.close(fig2)

