# Documentation for the AH-RJMCMC code 
## (Age-Hyperparameter Reverse-jump Markov Chain Monte Carlo)

**Version 1: March 27th 2018**  

Original release; produces the figures in the paper 
https://academic.oup.com/gji/article/215/3/2008/5101441.

**Version 2: Revised Jan 28th 2019**

Includes:  sampling of the ages from a normal distribution based on their current values, rather than directly from the prior. This speeds up convergence.
Increased flexibility in the type of data that can be read, including data files describing samples whose ages are a mixture of normal and uniform distributions; data types (in particular, brick samples with the same sample ID are not independent and are always defined with the same age); stratification of data.
Python version of the code is removed.

Authors: 
Phil Livermore (University of Leeds)
Alex Fournier (Institut de Physique du Globe de Paris, Paris)
Thomas Bodin (Universite Lyon, Lyon)

Maintained by Phil Livermore

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

`AH inputfile`

For example,

`AH input_Paris700`

runs the AH-RJMCMC code using the Paris700 dataset.

## Example datasets:

From the manuscript https://academic.oup.com/gji/article/215/3/2008/5101441

* Paris 700
* Hawaii
* Lubeck-Paris700

Additionally:

* Near East, Mixed type
* Near East, Group Level

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

to make a plot that shows the posterior ages of stratified data (specific to figure 16 of the referenced manuscript with the Lubeck-Paris dataset).

## Inputfile structure

Each inputfile contains all the information needed to run the AH-RJMCMC model.

Comment lines that begin with a '#' symbol are ignored.
Each line of information begins with a keyword (in either upper or lower case) followed by the associated information.
The information may appear in any order but there are no default values, so all information needs to be given.
The keywords that are needed to be specified in the input file are described below.

There are multiple input files in the Github directory that you can use as templates.


### Data_file: The relative path of the dataset to be read in 

### File_format (8 integer numbers) - These specify which columns of the data file should be read, in the order:  ID, age, delta age, intensity, delta intensity, data type, age distribution, stratification.

A data file may, for example, contain much more information than is required here and so only a subset of the columns need to be read.

ID : the column number of the sample ID (assumed to be text)

Age: the column number of the mean sample age

delta age: the column number for the age uncertainty; interpreted as the half-interval range for a uniformly distributed age or the standard deviation for a normally distributed age

intensity: the column number of the mean sample intensity

delta intensity: the column number of the standard deviation of the intensity

data type: the column number specifying the data type, or -1 if this not relevant. 
If -1 is specified, the data type is set internally within the code to 'O'(ther), which is never accessed.
In the data file, the choices of type are: P, B, C, S: pottery, brick, baked clay, slag.

age distribution: the column number specifying the age distribution flag or 
-1 if ALL the data ages are distributed uniformly (with centre point given by the age, and half-width given by the delta age)
-2 if ALL the data ages are distributed normally (with mean given by the age, and standard deviation given by delta age)

stratification: the column number specifying whether the data are stratified or 
-1 if stratification is to be ignored. 
In the data file, a value if 0 means that the datum is unstratified.
A value of i, where i>0, means the group of data with parameter i (there may be many such data) is assumed to represent a layer, within which the ages are unstratified. The next layer, of parameter i+1, is assumed to be another layer, but each datum within layer i must have an age less than (or greater than, depending on which way the ages are ordered) the data of parameter i+1.  If each layer contains only one datum, then these data are stratified and their ages must obey a strict ordering. A datum with a parameter of 0 delineates the boundary of the stratified data set. The age ordering of for the whole data is determined from the data with parameters 1 and 2. 

In all the above, the column index begins at 0.

A typical structure of File_format is:
File_format 0 2 3 4 5 -1 -1 -1  

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

sigma\_age: the standard deviation in age when proposing an alteration in the sample ages

### Age_frac: fraction of ages to change in a proposal 
(defined as the parameter beta in the manuscript referenced above).

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


