##### Functional Single Index Model - Compact Support Demonstration ######
##### Barinder Thind - 301193363 #####
######################################

### Loading libraries
library(devtools)
# install_github('YunlongNie/cFuSIM')
library(NonpModelCheck)
library(fda)
library(cFuSIM)

### Loading in the bike data
data(bike_cFuSIM)
bike

### Creating basis
timepts = bike$timepts
norder = 4 ## cubic B-spline
nbasis = norder + length(timepts) - 2; 
spline_basis = create.bspline.basis(rangeval=c(1,24),nbasis,norder,timepts)

### Turning data into functional data
wull = bike$temp
xfds = Data2fd(y = wull%>%t, argvals = bike$timepts)

### Pulling out the scalar response (total bike rentals over 102 weeks)
y = bike$y
train_sample = 1:length(y)
y = y[train_sample]
xfd = xfds[train_sample]

### Running the regression
res_c = cFuSIM_index(y, xfds, spline_basis, lambda = 10^5)
?cFuSIM_index

### Pulling out the beta estimate (index function)
beta_fd = fd(res_c$coefBeta, res_c$basisBeta)

### Plotting the index function 
fdagg(beta_fd,ylab="index function", xlab='time')

### Getting fitted values
score_fit = (res_c$score_fit)

### Getting predictions through the link function
pred_y = localpoly.reg(score_fit, y, degree.pol = 1, kernel.type = "gaussian",
                       bandwidth = "CV", deriv=0, points=score_fit)

?localpoly.reg

### plot the fitted integral vs the response 
plot(x = score_fit, y = y)
lines(pred_y$predicted[order(score_fit)], x=score_fit[order(score_fit)], col=4)

### compute MSE 
fity = pred_y$predicted
mse  = mean((y - fity)^2)

### Computing the fitted value
p = sum(res_c$coefBeta!=0)
bic = length(y)*log(mse)+log(length(y))*(p+1)


