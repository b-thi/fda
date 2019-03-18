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
spline_basis = create.bspline.basis(rangeval = c(1,24), nbasis = nbasis, 4, timepts)
bike_fd = Data2fd(timepts, t(bike$temp), spline_basis)

### Making list
bike_list = vector("list",2)
bike_list[[1]] = rep(1, 102)
bike_list[[2]] = bike_fd

### Creating a constant basis for the intercept
conbasis= create.constant.basis(c(1,24))

### Creating basis for beta
timepts = bike$timepts
spline_basis2 = create.bspline.basis(rangeval=c(1,24), 5, 4)

### List for coefficients
bike_beta_list  = vector("list",2)
bike_beta_list[[1]] = conbasis
bike_beta_list[[2]] = spline_basis2

### Fitting the model
lm_fda_bike = fRegress(y, bike_list, bike_beta_list)

bike_beta = lm_fda_bike$betaestlist[[2]]$fd

plot(bike_beta, xlab = "Time", ylab = "Average Bike Rentals")



### Getting predictons
bike_pred = lm_fda_bike$yhatfdobj
plot(bike_pred, y)
