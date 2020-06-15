# Script to make figures for manuscript
import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator,MaxNLocator)

def draw_composite(directory, symbols_file, title, ax1, ax2, credible, labels_font_size, plot_posterior_age_intensity=False, plot_spikes_times = False):
    
#0:  Tel hazor
#1:  Tel Megiddo
#2   Timna-30
#3   Judean jar handles
#4   Khirbat en-Nahas
#5   Qatna
#6   Other Syrian data
#7   Turkey Arslantepe
#8   Turkey, Kilise Tepe
#9   Turkey, Tell Tayinat

    #legend information
    color_edges = ["royalblue","royalblue","green","royalblue","royalblue","red","red","darkorange","darkorange","darkorange"]
    colors_faces = ["lightskyblue","lightskyblue","palegreen","lightskyblue","lightskyblue","tomato","tomato","gold","gold","gold"]
    err_colors = ["royalblue", "royalblue","green", "royalblue","royalblue","red","red","darkorange","darkorange","darkorange"]
    fmts = ["o","s", "o","^","X","o","^","o","s","^"]
    labels=["Israel, Tel Hazor","Israel, Tel Megiddo","Jordan, Khirbat en-Nahas","Israel, Timna-30","Israel, Judean jar handles","Syria, Qatna","Syria, other data","Turkey, Arslantepe","Turkey, Kilise Tepe","Turkey, Tell Tayinat"]

    
    
    lx, ly = np.loadtxt(directory +'credible_lower.dat', unpack=True)
    ux, uy = np.loadtxt(directory +'credible_upper.dat', unpack=True)
    median_x, median_y = np.loadtxt(directory +'median.dat', unpack=True)
    x, x_err, y, y_err, strat = np.loadtxt(directory + 'data.dat', unpack=True)
    
    ax1.fill_between(-lx, ly, uy, facecolor='plum', alpha=0.2, label='%i%% credible interval' % credible)
    ax1.plot(-lx,ly,'b--',linewidth=0.5)
    ax1.plot(-lx,uy,'b--',linewidth=0.5)
    
# load symbols data
    anom, atype, age, dt, Int, sd, vadm, sdvadm, val = np.loadtxt(symbols_file+".txt",unpack=True)
    anom=anom[(age<-500) & (age>-1200)]
    dt = dt[(age<-500) & (age>-1200)]
    Int=Int[(age<-500) & (age>-1200)]
    sd=sd[(age<-500) & (age>-1200)]
    val=val[(age<-500) & (age>-1200)]
    vadm=[(age<-500) & (age>-1200)]
    age=age[(age<-500) & (age>-1200)]
    
    if plot_posterior_age_intensity:
        new_age, new_intensity = np.zeros(len(age)), np.zeros(len(age))
        for index in range(len(age)):
            dist = np.loadtxt(directory+'Joint_distribution_data/Sample_'+str('{:04}'.format(1+index))+'.dat')
            age[index], Int[index], dt[index],sd[index] = np.median(dist[:,0]), np.median(dist[:,1]), np.std(dist[:,0]), np.std(dist[:,1])

    for i1 in range(0,10):
        if i1 == 2:  #reorder to show in order 0,1,2,3,4,... -> 0,1,3,4,2,5,6,7,8,9   #This is to get data from the same regions shown together
            i2 = 3
        elif i1==3:
            i2=4
        elif i1==4:
            i2 = 2
        else:
            i2 = i1
        
        mask=(val==i2)
        ax1.errorbar(-age[mask],Int[mask],yerr=sd[mask],xerr=dt[mask],markersize=6,fmt=fmts[i2], markeredgecolor=color_edges[i2],markerfacecolor=colors_faces[i2],ecolor=err_colors[i2],capsize=3,linewidth=0.5,alpha=1.0,label=labels[i2])

    ax1.plot(-median_x, median_y, 'navy', linewidth=2, label = 'Median', zorder=4)
    # Make x-tick labels invisible.
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.set_ylim(50,115)
    ax1.set_xlim(1250, 450)
    ax1.set_ylabel(r"Intensity$~(\mu$T)",fontfamily="serif",fontsize=labels_font_size,labelpad=8)
    ax1.yaxis.set_tick_params(labelsize=16)
    ax1.yaxis.set_minor_locator(MultipleLocator(5))
    
    # Make x-tick labels invisible.
    ax1.xaxis.set_ticks([])
    
    ax1.set_title(title, position=(0.01, 0.80),fontsize=14,loc='left',bbox=dict(facecolor=(1,1,1,1),edgecolor=(0,0,0,1)))
    ax1.spines['bottom'].set_visible(False)

    # dF/dt:
    plt.subplots_adjust(hspace=0.1, wspace=0)
    lx, ly = np.loadtxt(directory +'credible_dFdt_lower.dat', unpack=True)
    ux, uy = np.loadtxt(directory +'credible_dFdt_upper.dat', unpack=True)
    median_x, median_y = np.loadtxt(directory +'median.dat', unpack=True)
    
    ax2.plot(-lx,ly,'r--',linewidth=0.5)
    ax2.plot(-lx,uy,'r--',linewidth=0.5)
    ax2.fill_between(-lx, ly, uy, facecolor='coral', alpha=0.2, label='%i%% credible interval' % credible)
    dFdt = np.gradient(median_y,median_x)
    ax2.plot(-median_x, dFdt, 'red',linewidth=2, label = 'Median dF/dt')
    ax2.set_ylim(-1,1)
    
    duration_dFdt_threshold = 0.12
    dFdt_threshold = 0.60  # 5 x 0.12, and close to 0.62
    intensity_threshold = 5
    window_size = 5
    
    ax2.xaxis.set_tick_params(labelsize=16)
    ax2.yaxis.set_tick_params(labelsize=16)
    ax2.set_xlim(1250, 450)
    ax2.yaxis.set_ticks([-1,0,1])
    ax2.plot([1250,450],[dFdt_threshold,dFdt_threshold],'--k',linewidth=1)
    ax2.plot([1250,450],[-dFdt_threshold,-dFdt_threshold],'--k',linewidth=1)
    ax2.plot([1250,450],[duration_dFdt_threshold,duration_dFdt_threshold],':k',linewidth=1)
    ax2.plot([1250,450],[-duration_dFdt_threshold,-duration_dFdt_threshold],':k',linewidth=1)
    
    ax2.plot([1250,450],[0,0],'k',linewidth=1)
    ax2.set_ylabel(r'$dF/dt~(\mu T yr^{-1})$',fontfamily="serif",fontsize=labels_font_size,labelpad=14)
    
    ax2.yaxis.set_tick_params(labelsize=16)
    ax2.spines['top'].set_visible(False)
    ax2.yaxis.set_minor_locator(MultipleLocator(0.2))

    
    if plot_spikes_times:
        from spike_detector import find_spikes
        number_spikes, start_age, end_age, peak_age = find_spikes(median_x, median_y, window_size=window_size, intensity_threshold=intensity_threshold, dFdt_threshold=dFdt_threshold, duration_dFdt_threshold=duration_dFdt_threshold)
        for i in range(len(start_age)):
            ax1.axvspan(-start_age[i],-end_age[i],facecolor='lightskyblue',alpha=0.5,edgecolor='blue')
        
        #print('Maxima at age ',age[thispeak], ' has delta intensity of ', delta_intensity, ' and max dF/dt ',max_dFdt)
        print('number of spikes detected = ', number_spikes)

