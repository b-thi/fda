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
### for which to create the basis

# Fourier basis with 65 basis functions for a year
daybasis65 = create.fourier.basis(c(0,365), 65)

### Above, we need the number of basis functions and the entirety of period T
### 







###################### Functional Linear Model Code ###################### 

### Getting annual precipitation = Response
annualprec = log10(apply(daily$precav,2,sum))
ny = length(annualprec)

### Defining basis (Fourier)
tempbasis65  = create.fourier.basis(c(0,365),65)

### Smoothing the curve
## day.5 is the times for which the function is evaluated
## Then you pass the matrix of values
## Then you pass the temporary basis 
tempSmooth65 = smooth.basis(day.5, daily$tempav, tempbasis65)
tempfd65 = tempSmooth65$fd

### Now, need to create model matrix
## First is the vector of 1s
## Second is the smoothed matrix functional data object
templist      = vector("list",2)
templist[[1]] = rep(1,35)
templist[[2]] = tempfd65

### Need to create functional objects out of intercept and for beta
# create a constant basis for the intercept
conbasis   = create.constant.basis(c(0,365))

# Define the small number of basis functions for beta(t)
# We do not use roughness penalty here
nbasis = 11
betabasis5 = create.fourier.basis(c(0,365),nbasis)
betalist1  = vector("list",2)
betalist1[[1]] = conbasis
betalist1[[2]] = betabasis5

### Fitting the linear model
## First argument is the annual precipitation (reponse)
## Second argument is the functional data information
## Third argument is the betas
fRegressList1 = fRegress(annualprec, templist, betalist1)
names(fRegressList1)
betaestlist1  = fRegressList1$betaestlist
length(betaestlist1)

### Plotting the betas
# obtain beta(t)
tempbetafd1   = betaestlist1[[2]]$fd
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(tempbetafd1, xlab="Day", ylab="Beta for temperature")
### Due to uncontrollability of the beta, we need to add constraint so that
### it is smooth through the whole situation so we have to redo this all

### Doing predictions
## fitted value yhat
annualprechat1 = fRegressList1$yhatfdobj

## fitted residual
annualprecres1 = annualprec - annualprechat1

## sum squared residuals
SSE1 = sum(annualprecres1^2)

## sum squared residuals for the null model y = alpha + \epsilon
SSE0 = sum((annualprec - mean(annualprec))^2)

##### Doing penalized estimation regresion instead
### Essentially this means that we do it for the 
### betas when we do the basis functions for the beta
### Apparently, because this is periodic, the best
### choice here is likely the harmonic differential
### operator as opposed to some other operator

# Penalized Estimation

# Using the harmonic acceleration differential operator 
# to define roughness penalty on beta(t)
Lcoef = c(0,(2*pi/365)^2,0)
harmaccelLfd = vec2Lfd(Lcoef, c(0,365))

# We use 35 Fourier basis functions to represent beta(t)
##### Q - Why did we pick 35? instead of the 11 from before
betabasis35 = create.fourier.basis(c(0, 365), 35)

# Choosing Smoothing Parameters using cross-validation
loglam = seq(5,15,0.5)
nlam   = length(loglam)
SSE.CV = rep(NA,nlam)
for (ilam in 1:nlam) {
  print(paste("log lambda =", loglam[ilam]))
  lambda     = 10^(loglam[ilam])
  betalisti  = betalist1
  betalisti[[2]] = fdPar(betabasis35, harmaccelLfd, lambda)
  fRegi          = fRegress.CV(annualprec, templist, betalisti)
  SSE.CV[ilam]   = fRegi$SSE.CV
}

### Looking at plot of the cross-validated smoothness penalties
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(loglam, SSE.CV, type="b", lwd=2,
     xlab="log smoothing parameter lambda",
     ylab="Cross-validation score", cex.lab=2,cex.axis=2)

# Choose lambda which minimize SSE.CV
lambda      = 10^12.5
betafdPar.  = fdPar(betabasis35, harmaccelLfd, lambda)

betalist2      = betalist1
betalist2[[2]] = betafdPar.

# do functional linear model
annPrecTemp    = fRegress(annualprec, templist, betalist2)

# get beta(t)
betaestlist2   = annPrecTemp$betaestlist

# get fitted value yhat
annualprechat2 = annPrecTemp$yhatfdobj

# get the effective degrees of freedom
print(annPrecTemp$df)

#### Getting our beta back
# plot beta(t)
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(betaestlist2[[2]]$fd, xlab="Days",ylab="beta(t)",lwd=2,cex.lab=2,cex.axis=2)

# plot yhat vs. y: fitted vs. residuals
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(annualprechat2, annualprec, lwd=2,cex.lab=2,cex.axis=2)
abline(lm(annualprec ~ annualprechat2), lty='dashed', lwd=2)
abline(0,1,lty=1, lwd=2,col="red") ### What is this line?

##### Confidence Intervals

# fitted residuals
resid   = annualprec - annualprechat2

# estimate sigma^2
SigmaE. = sum(resid^2)/(35-annPrecTemp$df)
SigmaE  = SigmaE.*diag(rep(1,35))

# for smoothing temperature chat = a matrix * y
y2cMap  = tempSmooth65$y2cMap

# obtain point-wise standard error for beta(t)
stderrList = fRegress.stderr(annPrecTemp, y2cMap, SigmaE)

betafdPar      = betaestlist2[[2]]
betafd         = betafdPar$fd
betastderrList = stderrList$betastderrlist
betastderrfd   = betastderrList[[2]]

plot(betafd, xlab="Day", ylab="Temperature Reg. Coeff.",
     ylim=c(-6e-4,1.2e-03), lwd=2,cex.lab=2,cex.axis=2)
lines(betafd+2*betastderrfd, lty=2, lwd=2, col = "red")
lines(betafd-2*betastderrfd, lty=2, lwd=2, col = "red")
# The temperature has a significant positive effect on the annual precipitations from September to December. 
# There is no significant effect of the temperature on the annual precipitations in other time.


################# Redoing with some other data set (Height Data)

### Getting the response (male)
height_mean_m <- data.frame(growth$hgtm) %>% 
  summarise_all(funs(mean))

height_mean_m <- t(height_mean_m)[,1] + rnorm(39, mean = 5, sd = 5)


### Getting the response (female)
height_mean_f <- data.frame(growth$hgtf) %>% 
  summarise_all(funs(var))

height_mean_f <- t(height_mean_f)[,1]


### Creating fourier basis for males/females
height_basis  = create.bspline.basis(c(0,18), 4)

### Smoothing the curve
eval_years <- growth$age
height_smooth_m <- smooth.basis(eval_years, growth$hgtm, height_basis)
height_smooth_f <- smooth.basis(eval_years, growth$hgtf, height_basis)

### Pulling out the functional objects
height_fd_m = height_smooth_m$fd
height_fd_f = height_smooth_f$fd

### Now, need to create model matrix
## First is the vector of 1s
## Second is the smoothed matrix functional data object

# Initializing
h_m_list <- vector("list", 2)
h_f_list <- vector("list", 2)

# Intercepts
intercept_vec_m <- rep(1, ncol(growth$hgtm))
intercept_vec_f <- rep(1, ncol(growth$hgtf))

# Adding to list the functional objects
h_m_list[[1]] <- intercept_vec_m
h_f_list[[1]] <- intercept_vec_f
h_m_list[[2]] <- height_fd_m
h_f_list[[2]] <- height_fd_f

### Need to create functional objects out of intercept and for beta
# create a constant basis for the intercept
conbasis   = create.constant.basis(c(0,18))

# Define the small number of basis functions for beta(t)
# We do not use roughness penalty here
betabasis = create.bspline.basis(c(0,18), 4)
betalist1  = vector("list",2)
betalist1[[1]] = conbasis
betalist1[[2]] = betabasis

betabasis2 = create.bspline.basis(c(0,18), 4)
betalist2  = vector("list",2)
betalist2[[1]] = conbasis
betalist2[[2]] = betabasis2

### Fitting the linear model
## First argument is the annual precipitation (reponse)
## Second argument is the functional data information
## Third argument is the betas

# For females
fRegress_f <-  fRegress(height_mean_f, h_f_list, betalist1)
beta_est_f <- fRegressList1$betaestlist

### Plotting the betas
# obtain beta(t)
beta_fd_f <- beta_est_f[[2]]$fd
plot(beta_fd_f, xlab="Year", ylab="Beta for Time")

# For males
fRegress_m = fRegress(height_mean_m, h_m_list, betalist2)
beta_est_m  = fRegress_m$betaestlist

### Plotting the betas
# obtain beta(t)
beta_fd_m <-  beta_est_m[[2]]$fd
plot(beta_fd_m, xlab="Year", ylab="Beta for Time")


##### Doing penalized estimation regresion
### Essentially this means that we do it for the 
### betas when we do the basis functions for the beta
### Apparently, because this is periodic, the best
### choice here is likely the harmonic differential
### operator as opposed to some other operator

# Penalized Estimation

# Using the harmonic acceleration differential operator 
# to define roughness penalty on beta(t)
Lcoef = c(0,(2*pi/365)^2,0)
harmaccelLfd = vec2Lfd(Lcoef, c(0,18))

# We use 35 Fourier basis functions to represent beta(t)
##### Q - Why did we pick 35? instead of the 11 from before
betabasis35 = create.bspline.basis(c(0, 18), 4)

# Choosing Smoothing Parameters using cross-validation
loglam = seq(5,15,0.5)
nlam   = length(loglam)
SSE.CV = rep(NA,nlam)
for (ilam in 1:nlam) {
  print(paste("log lambda =", loglam[ilam]))
  lambda     = 10^(loglam[ilam])
  betalisti  = betalist1
  betalisti[[2]] = fdPar(betabasis35, harmaccelLfd, lambda)
  fRegi          = fRegress.CV(annualprec, templist, betalisti)
  SSE.CV[ilam]   = fRegi$SSE.CV
}

### Looking at plot of the cross-validated smoothness penalties
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(loglam, SSE.CV, type="b", lwd=2,
     xlab="log smoothing parameter lambda",
     ylab="Cross-validation score", cex.lab=2,cex.axis=2)

# Choose lambda which minimize SSE.CV
lambda      = 10^12.5
betafdPar.  = fdPar(betabasis35, harmaccelLfd, lambda)

betalist2      = betalist1
betalist2[[2]] = betafdPar.

# do functional linear model
annPrecTemp    = fRegress(annualprec, templist, betalist2)


xd
