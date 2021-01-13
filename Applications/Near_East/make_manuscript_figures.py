# Script to make figures for manuscript
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import os
import sys

if not os.path.exists('Outputs_Mixed_200M/input_file'):
    print('*'*50)
    print('Cannot find file: input_file')
    print('*'*50)
    sys.exit(0)


# Read some basic data from the input file
# This can be overwritten by either altering this file, or simply hardwiring the various parameters: e.g.
# age_min, age_max = 0, 100

for line in open('Outputs_Mixed_200M/input_file','r'):
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

# spikes definitions:
spikes_colour = ['red','orange','green','blue','magenta', 'black']
duration_dFdt_threshold = 0.12
dFdt_threshold = 0.6  # 5 x 0.12, and close to 0.62
intensity_threshold = 5
window_size = 5
from spike_detector import find_spikes

from matplotlib import colors
#print( colors.to_rgba('red') )
print('RGB of spike colours:')
for i in spikes_colour:
    print( colors.to_rgba(i) )

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator,MaxNLocator)
print('Building fig 1...')
files = ['Mixed_Plot', 'Fragment_Plot', 'Group-N2_Plot', 'Group-N3_Plot']
titles=[r"(a) Mixed level", r"(b) Thermal unit level",r"(c) Group level $N \geq 2$", r"(d) Group level, $N \geq 3$"]
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

# converts intensity to VADM (in units of 10^22 Am^2)
factor = 0.542

def int2VADM(x):
    return x/factor

def VADM2int(x):
    return x * factor

#legend information
color_edges = ["royalblue","royalblue","green","royalblue","royalblue","red","red","darkorange","darkorange","darkorange"]
colors_faces = ["lightskyblue","lightskyblue","palegreen","lightskyblue","lightskyblue","tomato","tomato","gold","gold","gold"]
err_colors = ["royalblue", "royalblue","green", "royalblue","royalblue","red","red","darkorange","darkorange","darkorange"]
fmts = ["o","s", "o","^","X","o","^","o","s","^"]
labels=["Israel, Tel Hazor","Israel, Tel Megiddo","Jordan, Khirbat en-Nahas","Israel, Timna-30","Israel, Judean jar handles","Syria, Qatna","Syria, other data","Turkey, Arslantepe","Turkey, Kilise Tepe","Turkey, Tell Tayinat"]

age_SHA, int_SHA, int_err_SHA = np.loadtxt('SHAWQ_IronAge_Palmyra.dat',unpack=True,usecols=(0,1,2),skiprows=1)
int_err_SHA = int_err_SHA * 2  #convert to 2 s.d.

for j in range(0,4):
    anom, atype, age, dt, Int, sd, vadm, sdvadm, val = np.loadtxt(files[j]+".txt",unpack=True)
    anom=anom[(age<-500) & (age>-1200)]
    dt = dt[(age<-500) & (age>-1200)]
    Int=Int[(age<-500) & (age>-1200)]
    sd=sd[(age<-500) & (age>-1200)]
    val=val[(age<-500) & (age>-1200)]
    vadm=[(age<-500) & (age>-1200)]
    age=age[(age<-500) & (age>-1200)]
    
    
    fig1,ax0 = plt.subplots(figsize=(9,3))
    for i1 in range(0,10):
        if i1 == 2:  #reorder to show in order 0,1,2,3,4,... -> 0,1,3,4,2,5,6,7,8,9
            i = 3
        elif i1==3:
            i=4
        elif i1==4:
            i = 2
        else:
            i = i1

        mask=(val==i)
        ax0.errorbar(-age[mask],Int[mask],yerr=sd[mask],xerr=dt[mask],markersize=6,fmt=fmts[i], markeredgecolor=color_edges[i],markerfacecolor=colors_faces[i],ecolor=err_colors[i],capsize=3,linewidth=0.5,alpha=1.0,label=labels[i])
    
    ax0.set_xlim(1210,490)
    ax0.set_ylim(45,115)

    ax0.xaxis.set_minor_locator(MultipleLocator(50))
    ax0.yaxis.set_minor_locator(MultipleLocator(5))

# Font et taille des légendes des axes
    ax0.set_xlabel("Age (BC)", fontfamily="serif",fontsize=16,labelpad=10)

    ax0.xaxis.set_tick_params(labelsize=16)
    ax0.yaxis.set_tick_params(labelsize=16)

# set up second axis for VADM:
    secay = ax0.secondary_yaxis('right', functions=(int2VADM, VADM2int))
    secay.yaxis.set_tick_params(labelsize=16)
    secay.yaxis.set_minor_locator(AutoMinorLocator())

    if j != 3:
        plt.setp(ax0.get_xticklabels(), visible=False)
        ax0.set_xlabel('')
    
    
    ax0.set_ylabel(r"Intensity$~(\mu$T)",fontfamily="serif",fontsize=16,labelpad=8)
    secay.set_ylabel(r'VADM$~(ZAm^2)$', fontfamily="serif",fontsize=16,labelpad=8)

        #ax0.set_ylabel(r"Archeointensity$~(\mu$T)",fontfamily="serif",fontsize=16,labelpad=8,color='white')
#secay.set_ylabel(r'VADM$~(ZAm^2)$', fontfamily="serif",fontsize=16,labelpad=8,color='white')

    secay.minorticks_on()

    if j ==0: # plot SCHA_DIF
        ax0.fill_between(-age_SHA, int_SHA-int_err_SHA, int_SHA+int_err_SHA, facecolor='lightgreen', alpha=0.2, edgecolor='green')
        ax0.plot(-age_SHA, int_SHA, 'darkgreen',lw=1)
        ax0.plot(-age_SHA, int_SHA-int_err_SHA,'g--',linewidth=0.5)
        ax0.plot(-age_SHA, int_SHA+int_err_SHA,'g--',linewidth=0.5)

    plt.title(titles[j], x=0.03, y = 0.85,fontsize=16,loc='left',bbox=dict(facecolor=(1,1,1,1),edgecolor=(0,0,0,1)))
    plt.text(0.82,0.09,r'$N_{data}$=%i'%len(age),fontsize=16,verticalalignment='center',transform=ax0.transAxes)
    plt.savefig(files[j]+'.pdf', bbox_inches='tight',pad_inches=0.0)

    plt.close()


# plot legend separately.

# remove errorbars:
from matplotlib import container
handles, labels = ax0.get_legend_handles_labels()
handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]

figLegend = plt.figure(figsize = (9,1))
plt.figlegend(handles,labels, loc = 'upper left',fontsize=12,labelspacing=0.4,ncol=3,frameon=False)
figLegend.savefig('legend.pdf')
plt.close()

# Make figure 2:
# Each part figure is created separately and compiled in LaTex. The figures are generated by a function, defined below.

print('Building fig 2...')
from draw_composite import draw_composite

dirs = ['Outputs_Frag_200M/', 'Outputs_Group-N2_200M/', 'Outputs_Group-N3_200M/']
titles=["(a) Thermal unit",r"(b) Group level, $N\geq 2$", r"(c) Group level, $N\geq 3$"]

Number_rows = 5
row_span_main_figure =3

for i in range(0,3):
    fig2 = plt.subplots(figsize=(10,4))
    ax1 = plt.subplot2grid((Number_rows, 1), (0,0), colspan=1,rowspan=row_span_main_figure)
    ax2 = plt.subplot2grid((Number_rows, 1), (row_span_main_figure,0), colspan=1,rowspan=Number_rows-row_span_main_figure)
    plt.subplots_adjust(hspace=0.1, wspace=0)
    draw_composite( dirs[i], files[i+1], titles[i], ax1, ax2, credible, 14, plot_posterior_age_intensity=True, plot_spikes_times = True )
   
   
    ax2.set_xlabel('Age (BC)',fontsize=14)
    
    if i == 0: #label the fragment plot with coloured stars
        median_x, median_y = np.loadtxt(dirs[0] +'median.dat', unpack=True)
        number_spikes_frag, start_age_frag, end_age_frag, peak_age_frag = find_spikes(median_x, median_y, window_size=window_size, intensity_threshold=intensity_threshold, dFdt_threshold=dFdt_threshold, duration_dFdt_threshold=duration_dFdt_threshold)
        for j in range(number_spikes_frag):
                ax1.plot(-peak_age_frag[j],108,'s',color='white', markersize=12)
                ax1.plot(-peak_age_frag[j],108,'*',color=spikes_colour[j], markersize=12)
        print('Figure 2: spikes found at ', peak_age_frag)
        print('Duration: ', np.array(end_age_frag) - np.array(start_age_frag))
    plt.savefig('Fig2_'+str(i)+'.pdf', bbox_inches='tight',pad_inches=0.0)
    plt.close()



print('Building fig 3...')
dir_base = "./"
dirs = [dir_base+'Outputs_Frag_min_3_200M/',dir_base+ 'Outputs_Frag25_min_3_200M/',dir_base+'Outputs_Frag_min_4_200M/',dir_base+ 'Outputs_Frag25_min_4_200M/',dir_base+ 'Outputs_Frag_min_5_200M/', dir_base+'Outputs_Frag25_min_5_200M/',dir_base+ 'Outputs_Frag_min_6_200M/', dir_base+'Outputs_Frag25_min_6_200M/']
titles=[r"(a) Minimum 3 $\mu$T",r"(b) Minimum 3 $\mu$T, 25 yr",r"(c) Minimum 4 $\mu$T",r"(d) Minimum 4 $\mu$T, 25 yr", r"(e) Minimum 5 $\mu$T",r"(f) Minimum 5 $\mu$T, 25 yr", r"(g) Minimum 6 $\mu$T",r"(h) Minimum 6 $\mu$T, 25 yr"]

Number_rows = 5
row_span_main_figure =3

for i in range(0,8):
    fig2 = plt.subplots(figsize=(10,5))
    ax1 = plt.subplot2grid((Number_rows, 1), (0,0), colspan=1,rowspan=row_span_main_figure)
    ax2 = plt.subplot2grid((Number_rows, 1), (row_span_main_figure,0), colspan=1,rowspan=Number_rows-row_span_main_figure)
    plt.subplots_adjust(hspace=0.1, wspace=0)
    draw_composite( dirs[i], files[1], titles[i], ax1, ax2, credible, 16, plot_posterior_age_intensity=True, plot_spikes_times = True)
    
    if i==6 or i==4:
        ax1.fill_between(-age_SHA, int_SHA-int_err_SHA, int_SHA+int_err_SHA, facecolor='lightgreen', alpha=0.2, edgecolor='green')
        ax1.plot(-age_SHA, int_SHA, 'darkgreen',lw=1)
        ax1.plot(-age_SHA, int_SHA-int_err_SHA,'g--',linewidth=0.5)
        ax1.plot(-age_SHA, int_SHA+int_err_SHA,'g--',linewidth=0.5)
    
    ax2.set_xlabel('Age (BC)',fontsize=14)
        
    plt.savefig('Fig3_'+str(i)+'.pdf', bbox_inches='tight',pad_inches=0.0)
    plt.close()

print('Building fig 4...')

dirs = [dir_base+'Outputs_Frag_no_Timna_age_model_Gaussian_200M/',dir_base+'Outputs_Frag_no_Timna_age_model_Gaussian_min_3_200M/',dir_base+'Outputs_Frag_no_Timna_age_model_Gaussian_min_5_200M/',dir_base+'Outputs_Frag_Extended_errors_uniform_no_Timna_age_model_Gaussian_200M/']
titles=[r"(a) Timna-30 ages Gaussian prior", r"(b) Timna-30 ages Gaussian prior, minimum intensity error 3 $\mu T$", r"(c) Timna-30 ages Gaussian prior, minimum intensity error 5 $\mu T$",r"(d) Timna-30 ages Gaussian prior, extended intensity errors"]

Number_rows = 5
row_span_main_figure =3

for i in range(0,len(titles)):
    fig2 = plt.subplots(figsize=(10,5))
    ax1 = plt.subplot2grid((Number_rows, 1), (0,0), colspan=1,rowspan=row_span_main_figure)
    ax2 = plt.subplot2grid((Number_rows, 1), (row_span_main_figure,0), colspan=1,rowspan=Number_rows-row_span_main_figure)
    plt.subplots_adjust(hspace=0.1, wspace=0)
    draw_composite( dirs[i], files[1], titles[i], ax1, ax2, credible, 16, plot_posterior_age_intensity=True, plot_spikes_times = True)
    
    ax1.fill_between(-age_SHA, int_SHA-int_err_SHA, int_SHA+int_err_SHA, facecolor='lightgreen', alpha=0.2, edgecolor='green')
    ax1.plot(-age_SHA, int_SHA, 'darkgreen',lw=1)
    ax1.plot(-age_SHA, int_SHA-int_err_SHA,'g--',linewidth=0.5)
    ax1.plot(-age_SHA, int_SHA+int_err_SHA,'g--',linewidth=0.5)

    ax2.set_xlabel('Age (BC)',fontsize=14)

    plt.savefig('Fig4_'+str(i)+'.pdf', bbox_inches='tight',pad_inches=0.0)
    plt.close()




print('Building fig 5...')


fig2,ax2 = plt.subplots(figsize=(6,4))
ax2.set_xlim(-60,60)
ax2.set_ylim(-20, 70)
ax2.set_xlabel('Centred age (yr)',fontsize=16)

median_x, median_y = np.loadtxt(dir_base+'Outputs_Frag_200M/' +'median.dat', unpack=True)
number_spikes, start_age, end_age, peak_age = find_spikes(median_x, median_y, window_size=window_size, intensity_threshold=intensity_threshold, dFdt_threshold=dFdt_threshold, duration_dFdt_threshold=duration_dFdt_threshold)
start_age, end_age = np.array(start_age), np.array(end_age)
print('Fragment level: found number of spikes :', number_spikes)

for j in range(number_spikes):
    # cut out spike between ages:
    spikes_profile = median_y[ (median_x >= start_age[j]) & (median_x <= end_age[j]) ]
    spikes_profile = spikes_profile - spikes_profile.max() - j * 10 + 50
    ax2.plot(np.arange(start_age[j],(end_age[j]+1))-peak_age[j], spikes_profile ,color=spikes_colour[j])
    #print( -np.arange(start_age[i],(end_age[i]+1)), median_y[ (median_x >= start_age[i]) & (median_x <= end_age[i]) ] )
    #plt.show()
ax2.set_title(r"Thermal unit spikes", position=(0.03, 0.85),fontsize=16,loc='left',bbox=dict(facecolor=(1,1,1,1),edgecolor=(0,0,0,1)))
ax2.xaxis.set_minor_locator(MultipleLocator(10))
#ax2.yaxis.set_minor_locator(MultipleLocator(10))
ax2.xaxis.set_tick_params(labelsize=16)
ax2.yaxis.set_tick_params(labelsize=16)
ax2.yaxis.set_ticks([])
ax2.spines['left'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.bar(-60,20,width=1,bottom = 0,color='k')
ax2.text(-60,20,r"$20~\mu$T",rotation='vertical',horizontalalignment='right',          verticalalignment='top',fontfamily="serif",fontsize=14)
fig2.savefig('Fig5_0.pdf', bbox_inches='tight',pad_inches=0.0)
plt.close()


# Plot shape of oldest spike as fn of intensity error budget
fig2,ax2 = plt.subplots(figsize=(6,4))
dirs2 = [dir_base+'Outputs_Frag_200M/',
        dir_base+'Outputs_Frag_min_3_200M/', dir_base+'Outputs_Frag_min_4_200M/', dir_base+'Outputs_Frag_min_5_200M/', dir_base+'Outputs_Frag_min_6_200M/']
titles=['Raw', r'Min 3$~\mu T$', r'Min 4$~\mu T$', r'Min 5$~\mu T$', r'Min 6$~\mu T$']
linestyle_fig2b = ["solid", "dashed" ,(0,(5,10)), "dotted", "dashdot",]
ax2.set_xlim(1020,950)
ax2.set_ylim(60, 110)
ax2.set_xlabel('Age BC (yr)',fontsize=16)
ax2.set_ylabel(r"Archeointensity$~(\mu$T)",fontfamily="serif",fontsize=16)

for j in range(len(titles)):
    median_x, median_y = np.loadtxt(dirs2[j] +'median.dat', unpack=True)
    number_spikes, start_age, end_age, peak_age = find_spikes(median_x, median_y, window_size=window_size, intensity_threshold=intensity_threshold, dFdt_threshold=dFdt_threshold, duration_dFdt_threshold=duration_dFdt_threshold)
    start_age, end_age = np.array(start_age), np.array(end_age)
    # age range: spike around 980 has index 1 for the raw dataset, and index 0 otherwise. With 6 micro T, there is no spike, so plot for the whole range.
    if j <= 2:
        start_for_plot =  start_age[1]
        end_for_plot = end_age[1]
    elif j == 4:
        start_for_plot =  -1020
        end_for_plot = -950
    else:
        start_for_plot =  start_age[0]
        end_for_plot = end_age[0]

#print( start_age, end_age, j)
    spikes_profile = median_y[ (median_x >= start_for_plot) & (median_x <= end_for_plot) ]
    ax2.plot(-np.arange(start_for_plot,(end_for_plot+1)), spikes_profile ,color=spikes_colour[1], linestyle=linestyle_fig2b[j],label=titles[j])
#print( -np.arange(start_age[i],(end_age[i]+1)), median_y[ (median_x >= start_age[i]) & (median_x <= end_age[i]) ] )
#plt.show()
#ax2.set_title(titles[i], position=(0.07, 0.85),fontsize=16,loc='left',bbox=dict(facecolor=(1,1,1,1),edgecolor=(0,0,0,1)))
#ax2.axvspan(-start_age_frag[spike_focus],-end_age_frag[spike_focus],facecolor='lightskyblue',alpha=0.5,edgecolor='blue')
fig2.legend(loc = 'upper center',fontsize=12,labelspacing=0.1, bbox_to_anchor=(0.32, 0.9),handlelength=5)
ax2.xaxis.set_tick_params(labelsize=16)
ax2.yaxis.set_tick_params(labelsize=16)


fig2.savefig('Fig5_1.pdf', bbox_inches='tight',pad_inches=0.0)
plt.close()


# Plot change in max dF/dt over spike duration
fig2,ax2 = plt.subplots(figsize=(2,7))
ax2.set_ylim(0,4)
number_cases = 6
ax2.set_xlim(-0.5,number_cases-0.5)
ax2.plot([-0.5,number_cases-0.5],[dFdt_threshold,dFdt_threshold],'--k',linewidth=1)

ax2.set_xticks(range(number_cases))
ax2.set_xticklabels(['Raw', 'Extended', r'Min 3$~\mu T$', r'Min 4$~\mu T$', r'Min 5$~\mu T$', r'Min 6$~\mu T$'], rotation='vertical')
# Pad margins so that markers don't get clipped by the axes
ax2.margins(0.2)
# Tweak spacing to prevent clipping of tick-labels
plt.subplots_adjust(bottom=0.15)

dirs = [dir_base+'Outputs_Frag_200M/', dir_base+'Outputs_Frag_extended_errors_normal_200M/',
         dir_base+'Outputs_Frag_min_3_200M/', dir_base+'Outputs_Frag_min_4_200M/', dir_base+'Outputs_Frag_min_5_200M/', dir_base+'Outputs_Frag_min_6_200M/']

median_x, median_y = np.loadtxt(dir_base+'Outputs_Frag_200M/' +'median.dat', unpack=True)
number_spikes, start_age, end_age, peak_age = find_spikes(median_x, median_y, window_size=window_size, intensity_threshold=intensity_threshold, dFdt_threshold=dFdt_threshold, duration_dFdt_threshold=duration_dFdt_threshold)
start_age, end_age = np.array(start_age), np.array(end_age)

for i in range(0,number_cases):
    median_x, median_y = np.loadtxt(dirs[i] +'median.dat', unpack=True)
    
    grad = np.gradient(median_y, median_x)
    if i % 2 == 1: ax2.axvspan(i-0.5,i+0.5,facecolor='silver',alpha=0.5)
    for myspike in range(number_spikes):
        offset = (myspike - number_spikes/2)/20
        spikes_max = max(np.abs(grad[ (median_x >= start_age[myspike]) & (median_x <= end_age[myspike]) ] ))
        ax2.plot(i+offset,spikes_max,'*',markersize=8,color=spikes_colour[myspike])
#print(i,myspike, spikes_max)
ax2.set_ylabel(r'Maximum $\vert dF/dt \vert$ during spike $(\mu T yr^{-1})$')
fig2.savefig('Fig5_2.pdf', bbox_inches='tight',pad_inches=0.0)
