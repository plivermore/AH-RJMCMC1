
import numpy as np
import sys
import os

"""
from numba import jit
@jit


def function1(age):
    n = len(age)
   # age = midpoint_age.copy()
    a =  np.array([ [ x * y for x in age ] for y in age ])
        
@jit

from numba import jit
@jit(nopython =True)
"""
def RJMCMC(Model_parameters, midpoint_age, delta_age, intensity, delta_intensity, stratification, Return_info ):
    """ Subroutine to compute the RJ-MCMC inversion for intensity

Inputs: 
    
    Model_parameters: a dictionary containing all the model parameters
    
    age, delta_age, intensity, delta_intensity, stratification: lists of data
    
    Return_info: a dictionary of output information
    
    """
    
# Seed the generator so we get the same values every time
    np.random.seed(seed  = 1) 
    
# set the best K to -10, in order that it is obvious if it is never updated.    
    k_best = -10
    
    k_max_array_bound = Model_parameters['K_max'] + 1;
#   Num_samples_to_store = int(np.ceil((Model_parameters['Nsamples']-Model_parameters['Burn_in'])/Model_parameters['thinning']))
    
    
#  Calculate number of collected samples for credible intervals -- if we are collecting.
    if Model_parameters['Calc_credible']:
        Num_samples_credible=int(np.ceil((Model_parameters['Nsamples']-Model_parameters['Burn_in'])*((100 - Model_parameters['credible'])/200.0)/Model_parameters['thinning']))  
        print('Collecting credible interval data' )
# Define an equally-spaced grid to define the model:
    X = np.linspace(Model_parameters['X_min'], Model_parameters['X_max'],Model_parameters['discretise_size'])

# predefine arrays to keep track of the credible intervals
    val_min, val_max = np.zeros(Model_parameters['discretise_size']), np.zeros(Model_parameters['discretise_size'])
    ind_min, ind_max = np.zeros(Model_parameters['discretise_size'],dtype=int), np.zeros(Model_parameters['discretise_size'],dtype=int)
    MINI, MAXI = np.zeros((Model_parameters['discretise_size'], Num_samples_credible)), np.zeros((Model_parameters['discretise_size'], Num_samples_credible))
    
# predefine other arrays   
    age, age_prop = np.zeros(len(midpoint_age)), np.zeros(len(midpoint_age))
    pt, pt_prop, pt_best  = np.zeros( (k_max_array_bound, 2)), np.zeros( (k_max_array_bound, 2)), np.zeros( (k_max_array_bound, 2))
    endpt = np.zeros(2)

# initialise working variables
    b = bb = AB = AD = PB = PD = ACV = PCV = AP = PP = PA = AA = 0

# Initialize - Define randomly the first model of the chain
    k = np.random.randint(Model_parameters['K_min'],high=Model_parameters['K_max']+1)

# set the data ages to be the given nominal age (i.e. discount any age error). 
# This is so datasets with stratification are valid for the initial model.
# If we randomised the ages, we'd have to check that stratification was satifisfied, 
# and it could take a while before we find a valid model.

    age = midpoint_age.copy()  #make a copy of the midpoint age.

# Check to ensure that the stratification constraints (if any) are satisifed
    if not check_stratification(age, stratification):
        print( 'INITIAL DATA SET IS NOT CONSISTENT WITH GIVEN STRATIFICATION CONSTRAINTS')
        sys.exit(0)

# Check to make sure that the ages do not extend past the model ends. For then we can't compute the likelihood.
# This only happens with normally distributed ages, for which the age can be any value with prob > 0.
#    age = np.array( [ max( Model_parameters['X_min'], min(a, Model_parameters['X_max'])) for a in age] )
    for i in range(len(age)):
        age[i] = max( Model_parameters['X_min'], min(age[i], Model_parameters['X_max']))

    for i in range(k):
        pt[i,0] = Model_parameters['X_min'] + np.random.rand() * (Model_parameters['X_max'] - Model_parameters['X_min']) #position of internal vertex
        pt[i,1] = Model_parameters['I_min'] + np.random.rand() * (Model_parameters['I_max'] - Model_parameters['I_min']) #magnitude of internal vertex
        
        endpt[0] = Model_parameters['I_min'] + np.random.rand() * (Model_parameters['I_max'] - Model_parameters['I_min'])
        endpt[1] = Model_parameters['I_min'] + np.random.rand() * (Model_parameters['I_max'] - Model_parameters['I_min'])
    
#  make sure the positions are sorted in ascending order based on age.
    #print(pt)
    #print('*')
    pt[0:k] = pt[pt[0:k,0].argsort()]
    #np.ndarray.sort(pt, axis = 0)
    #print(pt)
    
# COMPUTE INITIAL MISFIT
# suppress exp overflow warnings - this can happen at the early stages of the algorithm
    trash = np.seterr(over =  'ignore')
    
    like=0;
    interpolated_signal = Find_linear_interpolated_values( Model_parameters['X_min'], Model_parameters['X_max'], pt[0:k,:], endpt, age )
    #print(delta_intensity)
    
    #print( len(age))
    #print( intensity[81] )
    #print( interpolated_signal[81] )
    #print( delta_intensity[81] )
    #q = (intensity - interpolated_signal)**2 / (2.0 * delta_intensity**2)
    #print(q[81])
    if Model_parameters['running_mode'] == 1:
        like = np.sum( (intensity - interpolated_signal)**2 / (2.0 * delta_intensity**2) )
    else:
        like = 1.0
        

    like_best=like
    like_init=like
    print('Initial likelihood is %s' % like)

# setup output for model data
    if Model_parameters['output_model_freq'] > 0:
        output_models = open(os.path.join(Model_parameters['outputs_directory'],Model_parameters['output_model_name']), 'w')
        output_models.write('%d\n' % Model_parameters['discretise_size'])
        for i in range(Model_parameters['discretise_size']):
            output_models.write('%10.3f ' % X[i] )
        output_models.write('\n')

# setup output for joint distribution data
    if Model_parameters['output_joint_distribution_freq'] > 0:
        joint_distribution_directory = os.path.join(Model_parameters['outputs_directory'],'Joint_distribution_data')
        if not os.path.exists(joint_distribution_directory): os.makedirs(joint_distribution_directory)

        joint_dist = [0] * len(age)
        for i in range(len(age)):
            joint_dist[i] = open(os.path.join(joint_distribution_directory,'Sample_%04d.dat'% (i+1)),'w')
                 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%% START RJ-MCMC SAMPLING %%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for s in range(1,Model_parameters['Nsamples']+1):
        
# Print statistics of the chain.    
        if np.mod(s,Model_parameters['show'])==0 and s > Model_parameters['Burn_in']:
            print( 'Samples %d, Vertices %d, Acceptance: change F %7.2f, change age %7.2f, birth %7.2f, death %7.2f resample ages %7.2f, likelihood %8.2f' % (s,k,100.0*ACV/PCV if PCV != 0 else np.NaN,100.0*AP/PP if PP != 0 else np.NaN, 100.0*AB/PB if PB != 0 else np.NaN, 100.0*AD/PD if PD != 0 else np.NaN, 100.0*AA/PA if PA != 0 else np.NaN, like) )
         
        birth=move=death=change_age = change_value = 0
    
# initialise the proposed model with the current    
        age_prop = age.copy()
        pt_prop=pt.copy()
        endpt_prop = endpt.copy()
        like_prop = like
        k_prop = k
        prob = 1.0
        out = 1
    
#----------------------------------------------------------------------
# Every 3rd iteration, propose a new value
        if  np.mod(s,3)==0:
            if s>Model_parameters['Burn_in']:  PCV +=1
            change_value = 1
            k_prop = k
            ind = np.random.randint(0,high=k+2) #generate a random integer between 0 and k+1
            
# choose which interior point to change, and check bounds to see if outside prior

            if ind == k: # change left end point
                endpt_prop[0] = endpt[0] + np.random.randn() * Model_parameters['sigma_change']
                if endpt_prop[0] < Model_parameters['I_min'] or endpt_prop[0] > Model_parameters['I_max']:  out = 0
    
            elif ind == k+1: # change right end point
                endpt_prop[1] = endpt[1] + np.random.randn()*Model_parameters['sigma_change']
                if endpt_prop[1] < Model_parameters['I_min'] or endpt_prop[1] > Model_parameters['I_max']:  out = 0
    
            else: # change interior point
                #print(pt_prop[ind,1], pt[0,1])
                pt_prop[ind,1] += np.random.randn(1)*Model_parameters['sigma_change']
                #pt_prop[ind,1] = pt_prop[ind,1] + np.random.randn(1)*Model_parameters['sigma_change']
                #print(pt_prop[ind,1], pt[0,1])
                if pt_prop[ind,1] < Model_parameters['I_min'] or pt_prop[ind,1] > Model_parameters['I_max']: out = 0

# Every 3rd iteration iteration change the vertex positions
        elif np.mod(s,3)==1: # Change age position
            u = np.random.randint(0,high=3)  #choose randomly between 3 operations:
        
            if u == 0: # BIRTH ++++++++++++++++++++++++++++++++++++++
                birth=1
                if s> Model_parameters['Burn_in']: PB += 1
                k_prop = k+1
                #print(np.size(pt_prop), k_prop)
                pt_prop[k_prop-1,0]=Model_parameters['X_min'] + np.random.rand()*(Model_parameters['X_max']-Model_parameters['X_min'])
 # Ensure that the new age is different to all the others - if it is, set out = 0 and abandon this model
                if pt_prop[k_prop-1,0] in pt_prop[0:k_prop-1,0]: out = 0               
                if k_prop > Model_parameters['K_max']: out=0

# interpolate to find magnitude as inferred by current state
                if out == 1:
                    interpolated_signal = Find_linear_interpolated_values( Model_parameters['X_min'], 
                                                                          Model_parameters['X_max'], pt[0:k,:], endpt, pt_prop[k_prop-1,0] )
   
                    pt_prop[k_prop-1,1]=interpolated_signal+np.random.randn()*Model_parameters['sigma_birth']
            
# Get probability
                    prob=(1.0/(Model_parameters['sigma_birth']*np.sqrt( 2.0 * np.pi )) *
                    np.exp(-(interpolated_signal-pt_prop[k_prop-1,1])**2/(2.0*Model_parameters['sigma_birth']**2)) )
            
# Check BOUNDS to see if outside prior
            
                    if pt_prop[k_prop-1,1] > Model_parameters['I_max'] or  pt_prop[k_prop-1,1] < Model_parameters['I_min']:  out=0
                    if pt_prop[k_prop-1,0] > Model_parameters['X_max'] or  pt_prop[k_prop-1,0] < Model_parameters['X_min']:  out=0

# make sure the positions are sorted in ascending order.
                    pt_prop[0:k_prop] = pt_prop[pt_prop[0:k_prop,0].argsort()]
                    
            elif u == 1: # !  DEATH +++++++++++++++++++++++++++++++++++++++++
                death=1
                if s> Model_parameters['Burn_in']: PD += 1
    
                k_prop = k-1
                if k_prop < Model_parameters['K_min']: out=0
    
                if out == 1:
                    ind = np.random.randint(0,high=k) # choose a vertex between 0 and k-1
                    pt_death = pt[ind,:]
                    pt_prop = pt.copy()
                    pt_prop = np.delete(pt_prop,ind,axis=0) # remove point to be deleted
                    pt_prop = np.append( pt_prop, [[0,0]],axis=0)  #add row of zeros to end to make sure the shape doesn't change.
                    
# Get prob - interpolate 
                    interpolated_signal = Find_linear_interpolated_values( Model_parameters['X_min'], 
                                                                              Model_parameters['X_max'], pt_prop[0:k_prop,:], endpt_prop, pt_death[0] )
                    prob=( 1.0/(Model_parameters['sigma_birth']*np.sqrt(2.0*np.pi))  *  
                    np.exp(-(interpolated_signal -pt_death[1])**2/(2.0*Model_parameters['sigma_birth']**2))  )
                    

            else: # MOVE +++++++++++++++++++++++++++++++++++++++++++++++++++++++
                if s> Model_parameters['Burn_in']:  PP += 1
                move=1
                k_prop = k
                if k == 0: out = 0  #If there are no points to move, then we can't move any
                
                if out == 1: 
                    ind = np.random.randint(0,high=k) # choose a vertex between 0 and k-1
                    pt_prop[ind,0] = pt[ind,0]+np.random.randn()*Model_parameters['sigma_move']  #Normal distribution of move destination
                    if pt_prop[ind,0] < Model_parameters['X_min'] or pt_prop[ind,0] > Model_parameters['X_max']: out = 0 
                    
# Ensure that the new age is different to all the others - if it is, set out = 0 and abandon this model
                    if pt_prop[ind,0] in np.delete(pt[0:k],ind,axis=0): out = 0 


# make sure the positions are sorted in ascending order.
                    pt_prop[0:k_prop] = pt_prop[pt_prop[0:k_prop,0].argsort()]

        else: # every 3rd iteration change the ages
# select ages at random

            if s> Model_parameters['Burn_in']: PA += 1
            change_age = 1 
            num_age_changes = int(np.floor(len(age)/float(Model_parameters['age_frac'])))
            random_indices = np.random.randint(0,len(age),num_age_changes)
            for i in random_indices: #choose num_age_changes from the set of ages and perturb
                if Model_parameters['Age_distribution'] == 'U':
                    age_prop[i] = midpoint_age[i] + 2.0 * (np.random.rand(1)-0.5) * delta_age[i]
                else:
                    age_prop[i] = midpoint_age[i] + np.random.randn() * delta_age[i]
                if age_prop[i] < Model_parameters['X_min'] or age_prop[i] > Model_parameters['X_max']: out = 0
                

# Check to ensure that the stratification constraints (if any) are satisifed
                if not check_stratification(age_prop, stratification): out = 0
                    
# end: decide on what proposal to make

# COMPUTE MISFIT OF THE PROPOSED MODEL 

        if out==1:
            like_prop=0;
            interpolated_signal = Find_linear_interpolated_values( Model_parameters['X_min'], Model_parameters['X_max'], 
                                                                  pt_prop[0:k_prop,:], endpt_prop, age_prop )
            if Model_parameters['running_mode'] == 1:
                like_prop = np.sum( (intensity - interpolated_signal)**2 / (2.0 * delta_intensity**2) )
            else:
                like_prop = 1.0
        
       
    
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SEE WHETHER MODEL IS ACCEPTED
    
        accept=0
        alpha = 0
# The acceptance term takes different values according the the proposal that has been made.
        if out == 1:
            if birth==1:
    
                alpha = ((1.0/((Model_parameters['I_max']-Model_parameters['I_min'])*prob))*np.exp(-like_prop+like))
                if np.random.rand() <alpha:
                    accept=1
                    if s>Model_parameters['Burn_in']: AB += 1
            
            elif death==1:
                alpha = ((Model_parameters['I_max']-Model_parameters['I_min'])*prob)*np.exp(-like_prop+like)
                if np.random.rand() <alpha:
                    accept=1
                    if s>Model_parameters['Burn_in']: AD+=1
            
            else: # NO JUMP, i.e no change in dimension
                alpha = np.exp(-like_prop+like)
                if np.random.rand() <alpha:
                    accept=1
                    if s>Model_parameters['Burn_in']: 
                        if change_value == 1:
                            ACV += 1
                        elif move == 1:
                            AP += 1
                        elif change_age ==1:
                            AA += 1
                        else:
                            print('FATAL ERROR 1'); sys.exit(0)
                
#If accept, update the values
        if accept==1:
            k=k_prop
            pt=pt_prop.copy()
            like=like_prop
            endpt = endpt_prop.copy()
            age = age_prop.copy()
    

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Collect samples for the ensemble solution

        if s>Model_parameters['Burn_in'] and np.mod( s-Model_parameters['Burn_in'],Model_parameters['thinning'])==0:
            b+=1

       # if Model_parameters['joint_distribution_freq'] > 0:
        # write joint distribution data
    

    # CALL Find_linear_interpolated_values( k, x_min, x_max, pt, endpt, discretise_size, x(1:discretise_size), interpolated_signal)
  
    #IF( FREQ_WRITE_MODELS > 0) then
    #if( s>burn_in .AND. mod(s-burn_in,thin * FREQ_WRITE_MODELS) == 0) WRITE(15,format_descriptor) interpolated_signal(1:discretise_size)
    #ENDIF


# CALL Find_linear_interpolated_values( k, x_min, x_max, pt, endpt, discretise_size, x(1:discretise_size), interpolated_signal)
            interpolated_signal = Find_linear_interpolated_values( Model_parameters['X_min'], Model_parameters['X_max'], 
                                                              pt[0:k,:], endpt, X )
# DO THE AVERAGE
            Return_info['Av'] += interpolated_signal[:]

            
# build marginal distribution for ages:
            for i in range(len(age)):
                if Model_parameters['Age_distribution'] == 'U':
                    bin_index = int( np.floor( (age[i]-(midpoint_age[i]-delta_age[i])) / (delta_age[i] * 2.0) * Model_parameters['Nbins_age_marginal']))
                else:
                    bin_index = int( np.floor( (age[i]-(midpoint_age[i]-2.0 * delta_age[i])) / (delta_age[i] * 4.0) * Model_parameters['Nbins_age_marginal']))
# For normally distributed ages, bin centred on mean with a 2*standard deviation range each side.
# Should a value fall outside this range, then simply add to either the 1st or last bin.
                bin_index = max(bin_index, 0)
                bin_index = min(bin_index, Model_parameters['Nbins_age_marginal']-1)

                Return_info['Marginal_ages'][i,bin_index] += 1

# write model data to disk

            if Model_parameters['output_model_freq'] > 0:
                if np.mod( s-Model_parameters['Burn_in'], Model_parameters['thinning'] * Model_parameters['output_model_freq']) == 0:
                    for i in range(Model_parameters['discretise_size']):
                        output_models.write('%10.3f \n' % interpolated_signal[i] )

# collect joint distribution data
            if Model_parameters['output_joint_distribution_freq'] > 0 and np.mod( s-Model_parameters['Burn_in'], Model_parameters['thinning'] * Model_parameters['output_joint_distribution_freq']) == 0:
                interpolated_samples = Find_linear_interpolated_values( Model_parameters['X_min'], Model_parameters['X_max'], pt[0:k,:], endpt, age )
                for i in range(len(age)):
                    joint_dist[i].write('%15.3f %15.3f\n' % (age[i],interpolated_samples[i]) )
                        
# build marginal intensity density
            for i in range(len(X)):
                bin_index = int(np.floor( (interpolated_signal[i]-Model_parameters['I_min'])/ (Model_parameters['I_max']-Model_parameters['I_min']) * Model_parameters['Nbins']))
                if bin_index <0 or bin_index > Model_parameters['Nbins']-1:
                    print('FATAL ERROR, BIN_INDEX IS OUT OF RANGE')
                    print('MODEL POINT %s VALUE %s' %(i,interpolated_signal[i]) )
                    print('INTENSITY MIN/MAX %s %s ' %(Model_parameters['I_min'], Model_parameters['I_max'] ))
                    print('Model is %s %s %s' % (k,endpt,pt[0:k,:]))
                    print(age); print(''); print(interpolated_signal)
                    sys.exit(0)
                Return_info['Intensity_density'][i,bin_index] += 1

            
# Do (e.g.) the 95% credible interval by keeping the lowest and greatest 2.5% of
# all models at each sample point. We could either keep ALL the data and at
# the end determine these regions (but this is very costly in terms of memory), or keep a running list of the
# number of data points we need. At the end of the algorithm, simply take the
# maximum of the smallest points, and the min of the largest, to get the
# bounds on the credible intervals.
# Method:
# Num_samples_credible is the number of data points corresponding to 2.5% of the total number
# of samples (after thinning).
# Collect Num_samples_credible datapoints from the first Num_samples_credible samples.
# For each subsequent sample, see if the value should actually be inside
# the 2.5% tail. If it is, replace an existing value by the current value.
# Repeat.

            if Model_parameters['Calc_credible']:
                for i in range(Model_parameters['discretise_size']):
                    if b <= Num_samples_credible: 
                        #print(b-1, Num_samples_credible)
                        MINI[i,b-1]=interpolated_signal[i]
                        MAXI[i,b-1]=interpolated_signal[i]
                        if b == Num_samples_credible:
                            val_min[i] = MAXI[i,:].min(); ind_min[i] = MAXI[i,:].argmin()
                            val_max[i] = MINI[i,:].max(); ind_max[i] = MINI[i,:].argmax()

                    else:  #we've already filled the tails, now compare each data point to see whether it should be included or not.
                        if interpolated_signal[i] > val_min[i]:
                            MAXI[i,ind_min[i]] = interpolated_signal[i]
                            val_min[i] = MAXI[i,:].min(); ind_min[i] = MAXI[i,:].argmin()
                    
                        if interpolated_signal[i] < val_max[i]:
                            MINI[i,ind_max[i]] = interpolated_signal[i]
                            val_max[i] = MINI[i,:].max(); ind_max[i] = MINI[i,:].argmax()
                   
                    
# Build histogram of number of changepoints: k
            Return_info['Changepoint_hist'][k] += 1

# k can be zero here  - I think there is mistake in the fortran: k can never be zero.

# Do the histogram on change points
            for i in range(k):
                Return_info['Change_points'][bb]=pt[i,0]
                bb += 1

# process ALL models now...
        Return_info['Misfit'][s-1] = like
    
# Get the best model
        if like<like_best and accept == 1:
            pt_best = pt.copy()
            k_best = k
            endpt_best = endpt.copy()
            like_best = like
            age_best = age.copy()
    
# -----------------------------    
# end: the Sampling of the mcmc
# ----------------------------

    Return_info['Change_points'] = Return_info['Change_points'][0:bb] #only return non-zero values.
    Return_info['Av'] = Return_info['Av']/b
   # print(  Return_info['Intensity_density'][0,:], Return_info['Intensity_density'][10,:])

# Compute the credible intervals:
    Return_info['Credible_Sup'] = np.min ( MAXI[:,:], axis = 1)
    Return_info['Credible_Inf'] = np.max ( MINI[:,:], axis = 1)

# normalise marginal distributions
    Return_info['Intensity_density'][:,:] = np.array(Return_info['Intensity_density'][:,:]) / np.sum( Return_info['Intensity_density'][0,:] )
    
# Compute the mode
    Return_info['Mode'] = (0.5 + np.argmax(Return_info['Intensity_density'], axis=1)) / Model_parameters['Nbins'] *  (Model_parameters['I_max'] - Model_parameters['I_min']) + Model_parameters['I_min']

# Compute the median. Get the first instance of the count from the left being greater than half the total:
    for i in range(Model_parameters['discretise_size']):
        for j in range(Model_parameters['Nbins']):
            if np.sum ( Return_info['Intensity_density'][i,0:j]) >= np.sum( Return_info['Intensity_density'][i,:] )/2.0:
                #print(j, np.sum ( Return_info['Intensity_density'][i,0:j]), np.sum( Return_info['Intensity_density'][i,:] )/2.0); 
                Return_info['Median'][i] = (0.5 + j) / Model_parameters['Nbins'] *  (Model_parameters['I_max'] - Model_parameters['I_min']) + Model_parameters['I_min']
                break


# Calculate the "best" solution
    if k_best < 0:
        print('NO MINIMUM LIKELIHOOD SOLUTION FOUND')
        Return_info['Best'] = np.zeros(Model_parameters['discretise_size'])
    else:
        Return_info['Best'] = Find_linear_interpolated_values( Model_parameters['X_min'], Model_parameters['X_max'], 
                                                              pt_best[0:k_best,:], endpt_best, X )
# close file of model data
    if Model_parameters['output_model_freq'] > 0:
        output_models.close()
                    
# close file for joint distributions
    if Model_parameters['output_joint_distribution_freq'] > 0:
        for i in range(len(age)):
            joint_dist[i].close()
                        
    return


def check_stratification(age, stratification):
    """ Function check_stratification
    
        Returns 1 if all ages are consistent with stratification constraints
        Returns 0 if not.
    """
    age_check = [age[i] for i in range(len(age)) if stratification[i] == 1]
    sorted_ages = age_check.copy()
    sorted_ages.sort()
    
    #print(age_check)
    #print(sorted_ages)
    
    if sorted_ages == age_check: #ages are in increasing order
        return 1
    else:
        return 0
    


def Find_linear_interpolated_values(x_min, x_max, pt, endpt, query_ages):
    """
    Function: Find_linear_interpolated_values
        Find the linearly interpolated value using information from the end points
    
    """
#    import numpy as np
    k = np.size(pt,axis=0)
    linear_description_time, linear_description_intensity = np.zeros(k+2), np.zeros(k+2)
    linear_description_time[0], linear_description_intensity[0] = x_min, endpt[0] 
    linear_description_time[1:k+1], linear_description_intensity[1:k+1] = pt[:,0], pt[:,1]
    linear_description_time[k+1],linear_description_intensity[k+1]  = x_max, endpt[1]
    
    return np.interp(query_ages, linear_description_time, linear_description_intensity)
