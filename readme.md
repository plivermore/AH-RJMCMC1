# Documentation for the AH-RJMCMC code 
## (Age-Hyperparameter Reverse-jump Markov Chain Monte Carlo)

**Version 1: March 27th 2018**  

Original release; produces the figures in the paper 
https://academic.oup.com/gji/article/215/3/2008/5101441.

**Version 2: Revised July 2nd 2019**

Includes:  sampling of the ages from a normal distribution based on their current values, rather than directly from the prior. This speeds up convergence.
Increased flexibility in the type of data that can be read, including data files describing samples whose ages are a mixture of normal and uniform distributions; data types (in particular, brick samples with the same sample ID are not independent and are always defined with the same age); stratification of data.
Python version of the code is removed.
On-the-fly saving of some diagnostics to increase speed.

**Version 3: Revised July 24th 2019**
Includes: capability to invert also for a rescaled standard deviation of the intensity by an unknown constant (lambda), by including lambda in the list of hyper-parameters.

**Version 4: revised Jan 17th 2020**
Includes: an option for drawing from the prior age distributions when resampling the sample ages, and the ability to specify the number of age parameters to be resampled (rather than just a fraction of ages).

**Version 5: revised June 16th 2020**
Includes: an option for using a uniform distribution for the intensity, rather than normal distributions. 
Reproduces the outputs for the manuscript: "Archeomagnetic intensity variations during the era of geomagnetic spikes in the Near East".

**Version 6: revised March 19th 2021**
Includes: an option for linking stratified data sets


Authors: 
Phil Livermore (University of Leeds)
Alex Fournier (Institut de Physique du Globe de Paris, Paris)
Thomas Bodin (Universite Lyon, Lyon)

Maintained by Phil Livermore, Alex Fournier

Last update: March 19th 2021

## Overview

AH-RJMCMC is an open-source Fortran code that computes a transdimensional posterior probability distribution of archeomagnetic intensity as a function of time.

For speed, the main program is coded in Fortran 90; graphics are produced using Python (version 3).


## Compiling the code

To compile the code, type

make all

By default the Intel Fortran compiler ifort is used, but it is straightforward to edit line one of the Makefile to use a different compiler (the gfortran compiler option is given, but commented out). 

## Basic execution: running an example 

The code reads from an inputfile, which defines among other things, all the MCMC parameters, the directory where the outputs are to be written, and the location of the data file.
To run the Fortran code using a single inputfile, type

`AH Applications/AH-RJMCMC/inputfiles/inputfile`

For example, to run the method on the Paris-700 dataset type

`cd Applications/AH-RJMCMC`

`../../AH inputfiles/input_Paris700`


## Example datasets:

1) From the manuscript **https://academic.oup.com/gji/article/215/3/2008/5101441**

* Paris 700
* Hawaii
* Lubeck-Paris700

Additionally

* Paris-700, rescaled intensity error budget 

(See readme in the folder: Applications/AH-RJMCMC).

2) The study of spikes in the Near East: **https://doi.org/10.1016/j.pepi.2021.106657** 

* Near East, Mixed type
* Near East, Group Level, N=2,3
* Near East, Fragment level (original dataset, with modified variants)

(See the readme in the folder Applications/Near_East)

3) A linked sequence test case, see Applications/Linked_sequences 

## Making figures

From the outputs directory, run

`python ../make_plots.py`

for basic diagnostics. 
Run 

`python ../make_joint_distribution_plots.py 13`

to create a joint distribution plot of the age and intensity, for datum number 13.

Run

`python ../make_Posterior_with_shifted_ages.py`

to make a plot of the posterior mean data age and intensity overlaid with the original data. 

Run

`python ../make_Stratified_marginal_age_plot.py`

to make a plot that shows the posterior ages of stratified data (specific to figure 16 of the referenced GJI manuscript with the Lubeck-Paris dataset).

## Inputfile structure

Each inputfile contains all the information needed to run the AH-RJMCMC model.

Comment lines that begin with a '#' symbol are ignored.
Each line of information begins with a keyword (in either upper or lower case) followed by the associated information.
The information may appear in any order but there are no default values, so all information needs to be given.
The keywords that are needed to be specified in the input file are described below.

There are multiple input files in the Github directory that you can use as templates.


### Data_file: The relative path of the dataset to be read in 

### File_format (9 integer numbers): These specify which columns of the data file should be read, in the order below.  

* ID : the column number of the sample ID (assumed to be text)

* Age: the column number of the mean sample age

* delta age: the column number for the age uncertainty; interpreted as the half-interval range for a uniformly distributed age or the standard deviation for a normally distributed age

* intensity: the column number of the mean sample intensity

* delta intensity: the column number of the standard deviation of the intensity

* data type: the column number specifying the data type, or -1 if this not relevant. 
If -1 is specified, the data type is set internally within the code to 'O'(ther), which is never accessed.
In the data file, the choices of type are: P, B, C, S: pottery, brick, baked clay, slag.

* age distribution: 
	* the column number specifying the age distribution type. In the corresponding data file, the distribution should be specified by a column of single characters: "N" (normal) or "U" (uniform) for each datum. If uniform, the given age and delta age specifies the centre point and half-width; if normal, the given age and delta age specifies the mean and standard deviation.
	* -1 if ALL the data ages are distributed uniformly (with centre point given by the age, and half-width given by the delta age) or
	* -2 if ALL the data ages are distributed normally (with mean given by the age, and standard deviation given by delta age)

* stratification: the column number specifying stratification information or 
-1 if stratification is to be ignored. 
See the section below about denoting stratification within the data file.

* intensity likelihood: the column number specifying the intensity likelihood type or
-1 if all sample likelihoods normally distributed (this is the usual case).
If set determined by the data file, the likelihood should be defined as "N" or "U" in the corresponding row, for normal or uniform distribution respectively. The intensity error is taken to be, respectively, the standard deviation or the half interval in intensity.

In all the above, the column index begins at 0.

A typical structure of File\_format is:
File_format 0 2 3 4 5 -1 -1 -1 -1

### Burn_in: The number of samples for the burn-in period

### Nsamples: The number of samples overall (including burn-in)

### model_discretisation: the resolution in time of the returned posterior distribution
 
### Chain_parameters: show, thinning
show: the frequency of writing information to the screen
thinning: the frequency of using samples for diagnostics for the Markov chain
 

### running mode: 0 or 1. 
Switch between posterior and prior sampling. 
1 - normal, 
0 - set all likelihoods to 1 to recover the prior distributions.

### Age_bounds: Age Min and Max for model parametrisation

### Sigmas: sigma\_move, sigma\_change, sigma\_birth, sigma\_age 
Parameters that describe the various model perturbations.

sigma\_move: the standard deviation of a move in which the age of a vertex is altered.

sigma\_change: the standard deviation when a vertex is altered in intensity

sigma\_birth: the standard deviation in intensity when a new vertex is born

sigma\_age: the standard deviation in age when proposing an alteration in the sample ages.

If sigma_age is positive, interpreted as the standard deviation of the normal age perturbation (which is the same for all samples). If negative, the standard deviation of the normal age perturbation is interpreted as the fraction of the given age uncertainty; this differs between samples. If zero, ages are resampled according to prior distributions.
For example:

`Sigmas 30 2 2 3`: the ages resampled from a normal distribution with mean given by the current age and standard deviation of 3. 

`Sigmas 30 2 2 -0.5`: the ages resampled from a normal distribution with mean given by the current age and standard deviation of 0.5 * sd, where sd is the age uncertainty given in the supplied data file. 

`Sigmas 30 2 2 0`: the ages are resampled from the prior distribution 

### Age_frac: fraction of ages to change in a proposal 
(defined as the parameter beta in the manuscript referenced above).

### Num\_age\_changes: number of ages to change in a resample-age proposal 
Note that only one of age\_frac or num\_age\_changes should be set.

### Intensity\_prior: 
min/max bounds on the assumed uniformly distributed intensity prior for vertices: 
I\_min, I\_max in micro Tesla

### Num\_change\_points: K\_min, K\_max 
Number of internal vertices is bounded to be within the interval [K\_min, K\_max]

### Outputs\_directory: Directory for all outputs
This directory is created if it does not exist

### Credible: Credible interval 
e.g. 95 (for 95% credible interval) or any non-positive number for none

### Nbins: Number of bins for posterior marginals

### output\_model: Name and write frequency of model files. 
Enter a frequency of -1 for no models to be output.
For example, the line below causes the code to write every 10th model (after thinning) to the models file here "models.dat". 

`output_model models.dat 10`

### output\_joint\_distribution\_freq: Joint distribution output frequency
If frequency is -1 then nothing is output; otherwise this defines the frequency of the joint distribution after thinning.

### Optional parameters:
True_data:  relative path of "true" underlying evolution
If this is present, some of the plotting code will add this data to figures.

Plotting\_intensity\_range: min/max of plotting range of intensity
Used only by the plotting scripts. e.g.
`Plotting_intensity_range 50 80`

SEED: Seed the random number generator with seed integer n. The default value is 1.

Batch\_generate\_joint\_distributions: compile a batch of joint age/intensity plots using the given column of the datafile. In the datafile, 1=make plot, 0 = exclude. E.g.
`Batch_generate_joint_distributions 9`
uses column 9 (must be 0 or 1) of the datafile to generate multiple figures at once.
NB. This is entirely equivalent to generating each figure individually.

Sigma\_uncertain: used to indicate sampling parameters for the rescaled intensity errors.
Usage: `Sigma_uncertain uniform_bound sd_sigma sd_fraction`

uniform\_bound: the rescaling value (lambda) has a prior of U[0, uniform\_bound]

sd_sigma: at each proposal for which lambda is changed, the proposed value is taken from a normal distribution centred on the current value and whose standard deviation is sd\_sigma.

sd\_fraction: For sampling of the hyperparameters (sample ages and lambda) this indicates the relative frequency of lambda perturbations. Each time hyperparameters are resampled, lambda is perturbed with probability sd\_fraction (0 <= sd\_fraction <= 1) 

E.g.  `Sigma_uncertain 3 0.4 0.1`
assumes a prior of U[0,3] for lambda, uses perturbations of standard deviation 0.4, and the hyper-parameter change move alters lambda 10% of the time, and so the ages 90% of the time.

# The Data file
## Stratification
The column of stratification information within the data file can specify a variety of scenarios:

- no stratification (denoted 0)
- a single sequence of stratified layers (denoted 1,2,3,4...)
- multiple (independent) sequences of stratified layers (denoted 1a, 2a, 3a...1b,2b,3b,...etc)
- Sequences can be linked. For example, using 3a-1d, 4a, 5a....,1c,2c-1e,3c,4c.. (see below).

The integer and integer-letter format cannot be mixed: e.g. 1a, 2a, 3a, 4, 5, 6 is not allowed.


Multiple data may be present within a stratified layer. In such a case within each layer the data are assumed to be mutually independent but the layers themselves obey stratification.

The age ordering of for the whole data is determined from the data with parameters 1 and 2.

## Example 1:
The following data file describes a single sequence of stratified data

| Unused column 	| Unused column	| Age    	| Age error 	| Intensity 	| Intensity error 	| Stratification 	|
|---------------	|--------	|--------	|--------	|-----------	|-----------------	|----------------	|
| 0             	| 84     	| 1330.5 	| 47.5   	| 51.4      	| 0.7             	| 1              	|
| 0             	| 83     	| 1391.5 	| 108.5  	| 53.5      	| 1.0             	| 2              	|
| 0             	| 82     	| 1400.5 	| 99.0   	| 54.8      	| 1.9             	| 3              	|
| 0             	| 81     	| 1425   	| 124.0  	| 56.1      	| 1.4             	| 4              	|
| 0             	| 80     	| 1425   	| 124.0  	| 57.8      	| 2.9             	| 5              	|
| 0             	| 79     	| 1488.5 	| 60.5   	| 57.1      	| 2.0             	| 6              	|

 The data stratification is determined by data 1 and 2: age is increasing with stratification index.
 The ages are drawn according to the Monte-Carlo algorithm, provided that 
 
 **age datum 1 <= age datum 2 <= age datum 3 etc.**


## Example 2:
The following data file describes a multiple sequence of stratified data

| Unused column 	| Unused column 	| Age    	| Age error 	| Intensity 	| Intensity error 	| Stratification 	|
|---------------	|--------	|--------	|--------	|-----------	|-----------------	|----------------	|
| 0             	| 84     	| 1330.5 	| 47.5   	| 51.4      	| 0.7             	| 1a              	|
| 0             	| 83     	| 1391.5 	| 108.5  	| 53.5      	| 1.0             	| 2a              	|
| 0             	| 82     	| 1400.5 	| 99.0   	| 54.8      	| 1.9             	| 3a              	|
| 0             	| 82     	| 1212   	| 124.0  	| 56.1      	| 1.4             	| 1b              	|
| 0             	| 83     	| 1271   	| 124.0  	| 57.8      	| 2.9             	| 2b              	|
| 0             	| 43     	| 13034 	| 60.5   	| 57.1      	| 2.0             	| 3b              	|

 The data stratification is determined by the first data 1 and 2: age is increasing with stratification index.

The ages are drawn according to the Monte-Carlo algorithm, provided that 
 
 **age datum 1a <= age datum 2a <= age datum 3a**
 
 **age datum 1b <= age datum 2b <= age datum 3b**
 
## Example 3:
The following data file describes a single stratification sequence, each layer constraining multiple data

| Unused column 	| Unused 	| Age    	| Age error 	| Intensity 	| Intensity error 	| Stratification 	|
|---------------	|--------	|--------	|--------	|-----------	|-----------------	|----------------	|
| 0             	| 84     	| 1330.5 	| 47.5   	| 51.4      	| 0.7             	| 1              	|
| 0             	| 83     	| 1325 	| 108.5  	| 53.5      	| 1.0             	| 1              	|
| 0             	| 82     	| 1316 	| 99.0   	| 52.8      	| 1.9             	| 1              	|
| 0             	| 82     	| 1308   	| 124.0  	| 53.1      	| 1.4             	| 1              	|
| 0             	| 83     	| 1391   	| 124.0  	| 57.8      	| 2.9             	| 2              	|
| 0             	| 43     	| 1401 	| 60.5   	| 57.1      	| 2.0             	| 2              	|

 The data stratification is determined by the first data 1 and 2: age is increasing with stratification index.

The ages are drawn according to the Monte-Carlo algorithm, provided that 
 
 **all ages within layer 1 are <= all ages within layer 2**
 
 Within any layer there is no assumption of age ordering. Therefore, the first datum within layer 1 could have a drawn age greater or less than the second datum within layer 1. The only requirement in this example is that each datum within layer 1 has an age less than or equal to any datum within layer 2.

## Example 4: linked sequences
The following data file describes two linked sequence of stratified data

| Unused column 	| Unused column 	| Age    	| Age error 	| Intensity 	| Intensity error 	| Stratification 	|
|---------------	|--------	|--------	|--------	|-----------	|-----------------	|----------------	|
| 0             	| 84     	| 1330.5 	| 47.5   	| 51.4      	| 0.7             	| 1a              	|
| 0             	| 83     	| 1391.5 	| 108.5  	| 53.5      	| 1.0             	| 2a-1c              	|
| 0             	| 82     	| 1400.5 	| 99.0   	| 54.8      	| 1.9             	| 3a              	|
| 0             	| 82     	| 1395   	| 124.0  	| 56.1      	| 1.4             	| 1b              	|
| 0             	| 83     	| 1396   	| 124.0  	| 57.8      	| 2.9             	| 2b-3c              	|
| 0             	| 43     	| 1401 	| 60.5   	| 57.1      	| 2.0             	| 3b              	|
| 0             	| 43     	| 1393 	| 23.5   	| 59.1      	| 3.1             	| 1d-2c              	|
| 0             	| 43     	| 1420 	| 12.5   	| 62.1      	| 4.1             	| 2d             	|



The data stratification is determined by the first data 1 and 2: age is increasing with stratification index.

The ages are drawn according to the Monte-Carlo algorithm, provided that 
 
 **age datum 1a <= age datum 2a <= age datum 3a**
 
 **age datum 1b <= age datum 2b <= age datum 3b**
 
 **age datum 1d <= age datum 2d**
 
 Additionally, the three sequences (1a,2a,3a), (1b,2b,3b), (1d,2d) are linked by the condition:
 
 **age datum 1c <= age datum 2c <= age datum 3c**
 
 When specifying the linking, the indices do not need to be ordered (as in the example): here they come in the order (1c,3c,2c). 