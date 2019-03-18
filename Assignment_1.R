### Dr. Cao - Assignment 1 ###
### Barinder Thind ###
### Functional Data Analysis ###

### Loading Libraries
library(tidyverse)
library(fda)
library(refund)
library(fda.usc)
###

# loading data
data("melanoma")

# looking at data
str(melanoma)

## Preprocessing ##
## So, it looks like we will be estimating one functional
## observation here. The plan is to plot the data, decide
## on a basis, then decide on a smoothing method, then
## estimate the function. 
###################

# plotting data
melanoma %>% 
  ggplot( aes(x = year, y = incidence)) + 
  geom_point() + 
  theme_bw() +
  ggtitle("Melanoma Data") +
  xlab("Year") +
  ylab("Incidence Rate") +
  theme(plot.title = element_text(hjust = 0.5))

### a

# Choosing Smoothing Parameters using cross-validation
basis_val = seq(4, 35, 1)
len_basis = length(basis_val)
gcv_val_basis = c()

# Creating time vector
time_vec <- seq(1936, 1972, 1)

# Running loop
for (i in 1:len_basis) {
  
  # Particular value of lambda
  basis_choice = basis_val[i]
  
  # Picking nbasis
  bsplinebasis_test <- create.bspline.basis(rangeval = c(1936, 1972),
                                        norder = 4,
                                        nbasis = basis_choice)
  
  
  # Creating FD Object
  smooth_fd = smooth.basis(time_vec, melanoma$incidence,
                            bsplinebasis_test)
  
  # Pulling out error
  gcv_val_basis[i] = sum(smooth_fd$gcv)
  
}

# Plotting results
data.frame(basis = basis_val, gcv = gcv_val_basis) %>% 
  ggplot(aes(x = basis, y = gcv)) + 
  geom_line() +
  theme_bw() +
  xlab("Number of Basis") +
  ylab("Cross Validation Error") +
  ggtitle("Finding Appropriate Number of Basis") +
  geom_point(aes(x = basis[which.min(gcv)], 
                 y = gcv[which.min(gcv)]), 
             color = "red",
             size = 3) +
  theme(plot.title = element_text(hjust = 0.5))


### b

# Correct number of basis
bsplinebasis_final <- create.bspline.basis(rangeval = c(1936, 1972),
                                      norder = 4,
                                      nbasis = basis_val[which.min(gcv_val_basis)])

# evaluating basis to get matrix
evalsplinebasis_final <- eval.basis(time_vec, bsplinebasis_final)

# now let's figure out the coefficients using regression
chat <- (solve((t(evalsplinebasis_final))%*%evalsplinebasis_final))%*%(t(evalsplinebasis_final))%*%melanoma$incidence

# Doing by hand
t(chat)%*%t(evalsplinebasis_final)

str(chat)
str(t(evalsplinebasis_final))

plot(time_vec, t(chat)%*%t(evalsplinebasis_final), type = "l")


# now, we can create the functional object again
melFD = fd(chat, bsplinebasis_final,
           list("Year","", "Incidence"))
plot(melFD)

# plotting functional observation with original data
plotfit.fd(melanoma$incidence, melanoma$year, melFD, 
           lty=1, lwd=2, col=1)

# so, this might look better, but worried about overfitting
# hm, there seems to be an up and down effect, it could potentially,
# even make sense to use a fourier basis

### c

# Now, let's do the derivative version

# deriv 1
melFD_1 <- deriv.fd(melFD)
plot(melFD_1, 
     ylab = "First Derivative", 
     main = "Plotting the First Derivative")

# deriv 2
melFD_2 <- deriv.fd(melFD_1)
plot(melFD_2, 
     ylab = "Second Derivative",
     main = "Plotting the Second Derivative")

################ Question 2

### a

# creating bspline basis
bsplinebasis_q2 <- create.bspline.basis(rangeval = c(1936, 1972),
                                     norder = 4, nbasis = 12)

# now, using this basis, we evaluate it at the appropriate
# values of t
time_vec <- seq(1936, 1972, 1)

# Choosing Smoothing Parameters using cross-validation
lambda_val = seq(0, 10, 0.01)
len_lam = length(lambda_val)
gcv_val = c()

# Running loop
for (i in 1:len_lam) {
  
  # Particular value of lambda
  lambda = lambda_val[i]
  
  # Setting up
  fdParobj = fdPar(bsplinebasis_q2, 2, lambda)
  
  # Creating FD Object
  smoothlist = smooth.basis(time_vec, melanoma$incidence,
                            fdParobj)
  
  # Pulling out error
  gcv_val[i] = sum(smoothlist$gcv)
  
}

# Plotting results
data.frame(lambdas = lambda_val, gcv = gcv_val) %>% 
  ggplot(aes(x = lambdas, y = gcv)) + 
  geom_line() +
  theme_bw() +
  xlab("Lambda") +
  ylab("Cross Validation Error") +
  ggtitle("Finding Appropriate Lambda") +
  geom_point(aes(x = lambdas[which.min(gcv)], 
                 y = gcv[which.min(gcv)]), 
             color = "red",
             size = 3) +
  theme(plot.title = element_text(hjust = 0.5))

### b

# Creating best FD Object
melanoma_fd_rough = smooth.basisPar(time_vec, 
                                    melanoma$incidence, 
                                    bsplinebasis_q2, 
                                    lambda = lambda_val[which.min(gcv_val)])

# Plotting 
plotfit.fd(melanoma$incidence, melanoma$year, melanoma_fd_rough$fd, 
           lty=1, lwd=2, col=1)


### c

# deriv 1
melFD_rough_deriv <- deriv.fd(melanoma_fd_rough$fd)
plot(melFD_rough_deriv, 
     ylab = "First Derivative",
     main = "Plotting the First Derivative")

# deriv 2
melFD_rough_deriv2 <- deriv.fd(melFD_rough_deriv)
plot(melFD_rough_deriv2, 
     ylab = "Second Derivative",
     main = "Plotting the Second Derivative")
