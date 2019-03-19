### Dr. Cao - Assignment 3 ###
### Barinder Thind ###
### Functional Data Analysis ###

### Loading Libraries
library(tidyverse)
library(fda)
library(refund)
library(fda.usc)
###

# Loading data
load(file = "bike.RData")
bike$temp

bike$y

# Hourly temperature for 102 Saturdays
temp = t(bike$temp )
dim (temp)
timepts = bike$timepts
par(mfrow = c(1 ,1) , mar = c(8, 8, 4, 2))
matplot (timepts, temp, xlab = "Hours", ylab = "Temperature",
         cex.lab =1.5 , cex.axis =1.5 , type = "l")

# The total counts of casual bike rentals on Saturdays
rental = bike$y
length(rental)

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
lm_fda_bike = fRegress(rental, bike_list, bike_beta_list)

bike_beta = lm_fda_bike$betaestlist[[2]]$fd

plot(bike_beta, xlab = "Time", ylab = "Average Bike Rentals")



### Getting predictons
bike_pred = lm_fda_bike$yhatfdobj

### Plotting predictions
data.frame(bike_pred = bike_pred, y = y) %>% 
  ggplot(aes(x = bike_pred, y = y)) +
  geom_point() +
  theme_bw() +
  xlab("Bike Predictions (yhat)") + 
  ylab("Actual Values (y)") +
  ggtitle("Fitted vs. Actual") + 
  theme(plot.title = element_text(hjust = 0.5))
