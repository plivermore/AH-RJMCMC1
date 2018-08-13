import sys
import numpy as np
import RJMCMC_inversion
import os

# AH-RJMCMC code
# Written by Phil Livermore, Alex Fournier, Thomas Bodin
# March 27, 2018
#

if len(sys.argv) != 2:
    print("Syntax: python AH_RJMCMC.py <inputfile>")
    sys.exit(0)

# Read in inputfile - ignore all lines beginning with a "#" symbol

Model_parameters = {}

for line in open(sys.argv[1]):
    if not (line[0] == '#' or line == '\n'):
        if line.split()[0].upper() == 'Data_file'.upper():
            Data_filename = line.split()[1].rstrip(); 
        if line.split()[0].upper() == 'File_format'.upper():
            age_col, d_age_col, F_col, dF_col, Strat_col = int(line.split()[1]),int(line.split()[2]),int(line.split()[3]),int(line.split()[4]), int(line.split()[5])
    #
        if line.split()[0].upper() == 'Intensity_prior'.upper():
            Model_parameters['I_min'], Model_parameters['I_max'] =  float(line.split()[1]),float(line.split()[2])
        if line.split()[0].upper() == 'Burn_in'.upper():
            Model_parameters['Burn_in'] = int(line.split()[1])
        if line.split()[0].upper() == 'Nsamples'.upper():
            Model_parameters['Nsamples'] = int(line.split()[1])
        if line.split()[0].upper() == 'model_discretisation'.upper():
            Model_parameters['discretise_size'] = int(line.split()[1])
        if line.split()[0].upper() == 'Chain_parameters'.upper():
            Model_parameters['show'], Model_parameters['thinning'] = int(line.split()[1]), int(line.split()[2])
        if line.split()[0].upper() == 'Age_distribution'.upper():
            Model_parameters['Age_distribution'] = line.split()[1].rstrip().upper()
        if line.split()[0].upper() == 'running_mode'.upper():
            Model_parameters['running_mode'] = int(line.split()[1])
        if line.split()[0].upper() == 'Age_bounds'.upper():
            Model_parameters['X_min'], Model_parameters['X_max'] =  float(line.split()[1]),float(line.split()[2])
        if line.split()[0].upper() == 'Sigmas'.upper():
            Model_parameters['sigma_move'], Model_parameters['sigma_change'], Model_parameters['sigma_birth'] =  (float(line.split()[1]),float(line.split()[2]),
                                                                                                                  float(line.split()[3]) )
        if line.split()[0].upper() == 'Age_frac'.upper():
            Model_parameters['age_frac'] =  float(line.split()[1])
    #    if line.split()[0].upper() == 'Initial_change_points'.upper():
    #        Model_parameters['K_init'] =  int(line.split()[1])
        if line.split()[0].upper() == 'Num_change_points'.upper():
            Model_parameters['K_min'], Model_parameters['K_max'] =  int(line.split()[1]), int(line.split()[2])
        if line.split()[0].upper() == 'Nbins'.upper():
            Model_parameters['Nbins'] =  int(line.split()[1])
        if line.split()[0].upper() == 'output_model'.upper():
            Model_parameters['output_model_name'], Model_parameters['output_model_freq'] = line.split()[1], int(line.split()[2])
        if line.split()[0].upper() == 'output_joint_distribution_freq'.upper():
            Model_parameters['output_joint_distribution_freq'] = int(line.split()[1])
        if line.split()[0].upper() == 'Credible'.upper():
            Model_parameters['credible'] = float(line.split()[1])
        if line.split()[0].upper() == 'Outputs_directory'.upper():
            Model_parameters['outputs_directory'] = line.split()[1]


# Read data file
if Strat_col > 0:
    age, delta_age, intensity, delta_intensity, stratification = np.loadtxt(Data_filename,usecols=(age_col, d_age_col, F_col, dF_col, Strat_col), unpack=True, comments='#')
    stratification = [int(a) for a in stratification]
else:
    age, delta_age, intensity, delta_intensity = np.loadtxt(Data_filename,usecols=(age_col, d_age_col, F_col, dF_col), unpack=True, comments='#')
    stratification = np.zeros(len(age), dtype=int)

# check age range
if Model_parameters['Age_distribution'].upper() == 'U':
    if Model_parameters['X_min'] > min(age - delta_age) or Model_parameters['X_max'] < max(age+delta_age):
        print('Increase age range of model as it does not span all data points with uniform error bounds')
        print('Age range needs to include', min(age - delta_age),' - ', max(age + delta_age) )
        sys.exit(0)
else:
    if Model_parameters['X_min'] > min(age - 2.0 * delta_age) or Model_parameters['X_max'] < max(age+2.0 * delta_age):
        print('Increase age range of model as it does not span all data points with 2 sd in age')
        print('Age range needs to include', min(age - 2.0 * delta_age),' - ', max(age + 2.0 *  delta_age) )
        sys.exit(0)

if min( delta_intensity) == 0 or min( delta_age) == 0:
    print('Minimum error in either age or intensity is zero, incompatible with assumptions build into code')
    sys.exit(0)

if Model_parameters['K_max'] < Model_parameters['K_min']:
    print('Changepoints error: K_max >= K_min')
    sys.exit(0)

print( '-'*40)
print('Read ', len(age),' data items')
print('Data age range is ', min(age), ' TO ', max(age))
print('Model age range is ', Model_parameters['X_min'], ' TO ', Model_parameters['X_max'])
print('Running mode is: ', 'Posterior sampling' if Model_parameters['running_mode'] == 1 else 'Prior sampling')
print('No stratification constraints active' if max( stratification ) == 0 else 'Stratification constraints active' )
print('Outputs directory is ' + Model_parameters['outputs_directory'])
print( '-'*40)

# Define an equally-spaced grid to define the model:
X = np.linspace(Model_parameters['X_min'], Model_parameters['X_max'],Model_parameters['discretise_size'])
if Model_parameters['credible'] > 0: 
    Model_parameters['Calc_credible'] = True
else:
    Model_parameters['Calc_credible'] = False
    
# Run the RJMCMC model
Return_info = {}
# Setup output structures

Num_samples_to_store = int(np.ceil((Model_parameters['Nsamples']-Model_parameters['Burn_in'])/Model_parameters['thinning']))
k_max_array_bound = Model_parameters['K_max'] + 1

Model_parameters['Nbins_age_marginal'] = Model_parameters['Nbins']

Return_info['Av'] = np.zeros(Model_parameters['discretise_size'])
Return_info['Mode'] = np.zeros(Model_parameters['discretise_size'])
Return_info['Median'] = np.zeros(Model_parameters['discretise_size'])
Return_info['Credible_Sup'] = np.zeros(Model_parameters['discretise_size'])
Return_info['Credible_Inf'] = np.zeros(Model_parameters['discretise_size'])
Return_info['Best'] = np.zeros(Model_parameters['discretise_size'])
Return_info['Change_points'] = np.zeros( Num_samples_to_store * k_max_array_bound )
Return_info['Misfit'] = np.zeros( Model_parameters['Nsamples'])
Return_info['Marginal_ages'] = np.zeros((len(age),Model_parameters['Nbins_age_marginal'] ),dtype=int )
Return_info['Changepoint_hist'] = np.zeros( k_max_array_bound,dtype=int )
Return_info['Intensity_density'] = np.zeros((Model_parameters['discretise_size'],Model_parameters['Nbins'] ) )



if not os.path.exists(Model_parameters['outputs_directory']):
    os.makedirs(Model_parameters['outputs_directory'])

RJMCMC_inversion.RJMCMC(Model_parameters, age, delta_age, intensity, delta_intensity, stratification, Return_info)

print('Writing outputs to disk...')
np.savetxt(os.path.join(Model_parameters['outputs_directory'],'data.dat'), np.column_stack((age, delta_age, intensity, delta_intensity,stratification)),fmt = '%10.3f')
np.savetxt(os.path.join(Model_parameters['outputs_directory'],'k_histogram.dat'), np.column_stack((range(0,Model_parameters['K_max']+1),Return_info['Changepoint_hist'])),fmt = '%i %i')
np.savetxt(os.path.join(Model_parameters['outputs_directory'],'credible_upper.dat'), np.column_stack((X,Return_info['Credible_Sup'])),fmt = '%10.3f')
np.savetxt(os.path.join(Model_parameters['outputs_directory'],'credible_lower.dat'), np.column_stack((X,Return_info['Credible_Inf'])),fmt = '%10.3f')
np.savetxt(os.path.join(Model_parameters['outputs_directory'],'average.dat'), np.column_stack((X,Return_info['Av'])),fmt = '%10.3f')
np.savetxt(os.path.join(Model_parameters['outputs_directory'],'best_fit.dat'), np.column_stack((X,Return_info['Best'])),fmt = '%10.3f')
np.savetxt(os.path.join(Model_parameters['outputs_directory'],'mode.dat'), np.column_stack((X,Return_info['Mode'])),fmt = '%10.3f')
np.savetxt(os.path.join(Model_parameters['outputs_directory'],'median.dat'), np.column_stack((X,Return_info['Median'])),fmt = '%10.3f')
np.savetxt(os.path.join(Model_parameters['outputs_directory'],'misfit.dat'), np.column_stack((range(1,Model_parameters['Nsamples'],100), Return_info['Misfit'][::100])),fmt = '%i %10.3f')
np.savetxt(os.path.join(Model_parameters['outputs_directory'],'changepoints.dat'), Return_info['Change_points'],fmt = '%10.3f')

# save the input file to the outputs directory
from shutil import copyfile
file = os.path.join(Model_parameters['outputs_directory'],'input_file')
copyfile(sys.argv[1], file)

# write intensity density
f = open(os.path.join(Model_parameters['outputs_directory'],'intensity_density.dat'),'w')
f.writelines("%s %s\n" % (Model_parameters['discretise_size'], Model_parameters['Nbins']))
for i in range(Model_parameters['discretise_size']):
    for j in range(Model_parameters['Nbins']):
        f.writelines("%10.3f %10.5f %10.3f\n" % (X[i] , (j-0.5)/Model_parameters['Nbins'] * (Model_parameters['I_max']-Model_parameters['I_min']) + Model_parameters['I_min'],
        Return_info['Intensity_density'][i,j] ))
f.close()


# write age marginals
f = open(os.path.join(Model_parameters['outputs_directory'],'age_marginals.dat'),'w')
f.writelines("%s %s\n" % (len(age), Model_parameters['Nbins_age_marginal']))
for i in range(len(age)):
    for j in range(Model_parameters['Nbins_age_marginal']):
        if Model_parameters['Age_distribution'] == 'U':
            bin_age = age[i]-delta_age[i] + float(j-1)/Model_parameters['Nbins_age_marginal'] * 2 * delta_age[i]
        else:
            bin_age = age[i]-2 * delta_age[i] + float(j-1)/Model_parameters['Nbins_age_marginal'] * 4 * delta_age[i]
        
        f.writelines("%10.3f %12.5f\n" % ( bin_age, Return_info['Marginal_ages'][i,j] ))
f.close()
