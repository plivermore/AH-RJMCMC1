# Script to combine the mean and credible intervals into a single file
import numpy as np

time, av = np.loadtxt('average.dat', unpack=True)
lx, ly = np.loadtxt('credible_lower.dat', unpack=True)
ux, uy = np.loadtxt('credible_upper.dat', unpack=True)

f= open("combined_mean_credible.txt","w+")
for index in range(len(time)):
	f.write("%f %f %f %f\n" %(time[index],av[index], ly[index], uy[index]))
f.close()


