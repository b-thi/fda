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

# plotting
ggplot(data = melanoma, aes(x = year, y = incidence)) + geom_point() + theme_bw()
## Looks like a bspline basis might be reasonable here

# creating bspline basis
bsplinebasis <- create.bspline.basis(rangeval = c(1936, 1972),
                                 norder = 2)

# now, using this basis, we evaluate it at the appropriate
# values of t
time_vec <- seq(1936, 1972, 1)

# evaluating basis
evalsplinebasis <- eval.basis(time_vec, bsplinebasis)

# now let's figure out the coefficients
chat <- (solve((t(evalsplinebasis))%*%evalsplinebasis))%*%(t(evalsplinebasis))%*%melanoma$incidence

# now, we can create the functional object
melFD = fd(chat, bsplinebasis,
            list("Year","", "Incidence"))

# plotting
plotfit.fd(melanoma$incidence, melanoma$year, melFD, 
           lty=1, lwd=2, col=1)

# Looking at the above plot, it actually looks like
# there could be reason to use a cubic spline basis,
# with a few basis functions, let's try this instead

# creating bspline basis (order 4)
bsplinebasis2 <- create.bspline.basis(rangeval = c(1936, 1972),
                                     norder = 4,
                                     nbasis = 15)

# evaluating basis
evalsplinebasis2 <- eval.basis(time_vec, bsplinebasis2)

# now let's figure out the coefficients using regression
chat2 <- (solve((t(evalsplinebasis2))%*%evalsplinebasis2))%*%(t(evalsplinebasis2))%*%melanoma$incidence


# now, we can create the functional object again
melFD2 = fd(chat2, bsplinebasis2,
           list("Year","", "Incidence"))

# plotting again
plotfit.fd(melanoma$incidence, melanoma$year, melFD2, 
           lty=1, lwd=2, col=1)

# so, this might look better, but worried about overfitting
# hm, there seems to be an up and down effect, it could potentially,
# even make sense to use a fourier basis

# Now, let's do the derivative version

# deriv 1
melFD2_1 <- deriv.fd(melFD2)
plot(melFD2_1)

# deriv 2
melFD2_2 <- deriv.fd(melFD2_1)
plot(melFD2_2)

################ Question 2

# creating bspline basis
bsplinebasis_q2 <- create.bspline.basis(rangeval = c(1936, 1972),
                                     norder = 4, nbasis = 15)

# now, using this basis, we evaluate it at the appropriate
# values of t
time_vec <- seq(1936, 1972, 1)

# Defining differential operator
Lcoef = c(0,(2*pi/365)^2,0)
harmaccelLfd = vec2Lfd(Lcoef, c(1936, 1972))

# Creating variable with parameter
melanoma_par <- fdPar(bsplinebasis_q2, 2, 1)

# Smoothing and pulling out functional obs
melFD_smooth <- smooth.basis(time_vec, melanoma$incidence,
                        melanoma_par)$fd

# Plotting
plotfit.fd(melanoma$incidence, melanoma$year, melFD_smooth, 
           lty=1, lwd=2, col=1)




# Choosing Smoothing Parameters using cross-validation
loglam = seq(0, 10, 0.01)
nlam = length(loglam)
dfsave = rep(NA,nlam)
gcvsave = rep(NA,nlam)

# Running loop
for (ilam in 1:nlam) {
  cat(paste("log10 lambda =",loglam[ilam],"\n"))
  lambda = loglam[ilam]
  fdParobj = fdPar(bsplinebasis_q2, 2, lambda)
  smoothlist = smooth.basis(time_vec, melanoma$incidence,
                            fdParobj)
  dfsave[ilam] = smoothlist$df
  gcvsave[ilam] = sum(smoothlist$gcv)
}

plot(loglam, gcvsave)


test = smooth.basisPar(time_vec, melanoma$incidence, bsplinebasis_q2, lambda = loglam[which.min(gcvsave)])

plotfit.fd(melanoma$incidence, melanoma$year, test$fd, 
           lty=1, lwd=2, col=1)

