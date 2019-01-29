###### Practice FDA ######

### Libraries
library(fda)

### Defining a fourier basis
# Below, we have the end points of observations
# And the number of basis functions we are using to
# define the functional observation
daybasis65 = create.fourier.basis(c(0,365), 65)

### Defining a b-spline basis
# Below, degree = 12 and over 9 equally spaced
# points from 0 to 10
splinebasis = create.bspline.basis(c(0,10), 13)
plot(splinebasis)

### Creating basis for height data
# Has good second derivatives (order = 6)
# Defined with knots at the ages
heightbasis = create.bspline.basis(c(1, 18), 35, 6, growth$age)
summary(heightbasis)

### Evaluating basis at a set of values of t
# Not super sure what this is yet but I know
# what this looks like - possibly used in coefficient
# calculation? I think that this is like we know when
# he functions are evaluated so we need to know those values
# for the purpose of the linear combination - we will know
# the value I think at some time t and just multiply
# any row of the below matrix by the coefficients atfer
# we get the coefficients by calculating them from the
# actual data points
tvec = c(1:18)
eval.basis(tvec, heightbasis)


