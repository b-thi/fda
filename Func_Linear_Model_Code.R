##### Functional Linear Model - Demonstration ######
##### Barinder Thind - 301193363 #####
######################################

### Loading libraries
library(fda)

### Getting Scalar Response
y = bike$y
n = length(y)

### Creating basis
timepts = bike$timepts
norder = 4 
nbasis = norder+length(timepts)-2; 
spline_basis = create.bspline.basis(rangeval=c(1,24), nbasis, norder, timepts)
timepts = timepts + rnorm(24, mean = 0, sd = 0.1)
timepts[1] = 1
timepts[24] = 24
bike_fd = Data2fd(timepts, t(bike$temp), spline_basis)

### Making list
bike_list = vector("list",2)
bike_list[[1]] = rep(1, 102)
bike_list[[2]] = bike_fd

### create a constant basis for the intercept
conbasis   = create.constant.basis(c(1,24))

### Creating basis for beta
timepts = bike$timepts
norder = 4 
nbasis = norder + length(timepts) - 2; 
spline_basis2 = create.bspline.basis(rangeval=c(1,24), nbasis, norder, timepts)

### List for coefficients
bike_beta_list  = vector("list",2)
bike_beta_list[[1]] = conbasis
bike_beta_list[[2]] = spline_basis2

### Fitting the model
lm_fda_bike = fRegress(y, bike_list, bike_beta_list)

bike_list


# fit the functional linear model 
fRegressList1 = fRegress(annualprec,templist,betalist1)
names(fRegressList1)
betaestlist1  = fRegressList1$betaestlist
length(betaestlist1)
# betaestlist1 has two elements. The first element is the intercept
# The second element is the slope beta(t)

# obtain beta(t)
tempbetafd1   = betaestlist1[[2]]$fd

quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(tempbetafd1, xlab="Day", ylab="Beta for temperature")
