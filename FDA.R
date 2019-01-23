####################################
###### FDA From the Beginning ######
###### Barinder Thind - 301193363 ##
####################################

####################################
####################################
####################################

### Loading libraries
library(fda)
library(fda.usc)
library(tidyverse)

### First, we need to create a basis for the functional observations
### that we want to make, then we need to evaluate that basis at some particular
### time points. Also, note that we need to pick a number of basis functions
### for which to create the 

# Fourier basis with 65 basis functions for a year
daybasis65 = create.fourier.basis(c(0,365), 65)
