##### SCAD Function #####

### Defining Function
p_scad <- function(beta, lambda, a){
  
  vals = c()
  
  for (i in 1:length(beta)) {
    
    if (0 <= beta[i] & beta[i] <= lambda){
      vals[i] <- lambda*beta[i]
    }
    
    if (lambda < beta[i] & beta[i] < a*lambda){
      vals[i] <- -(beta[i]^2 - 2*a*lambda*beta[i] + lambda^2)/(2*(a - 1))
    }
    
    if (beta[i] >= a*lambda){
      vals[i] <- ((a + 1)*lambda^2)/2
    }
    
  }

  return(as.numeric(vals))
  
}

### Defining beta values
beta_func <- function(nums){
  vals <- nums^2
  return(vals)
}

nums = seq(0, 10, 0.001)

### Running function
p_scad_val <- p_scad(beta_func(nums), 0.1, 3.7)

plot(beta_func(nums), p_scad_val)


### Integrating
integral_vals = c()
for (i in 1:100) {
 integral_vals[i] <- integrate(p_scad, a = 3.7, lambda = 70, 0, i)[[1]]/i
}

plot(1:100, integral_vals, type = "l")

integral_vals = c()
for (i in 1:100) {
  integral_vals[i] <- integrate(p_scad, a = 3.7, lambda = 3, 0, i)[[1]]/i
}

lines(1:100, integral_vals, col = "red")

integral_vals = c()
for (i in 1:100) {
  integral_vals[i] <- integrate(p_scad, a = 3.7, lambda = 0.5, 0, i)[[1]]/i
}

lines(1:100, integral_vals, col = "blue")

integral_vals = c()
for (i in 1:100) {
  integral_vals[i] <- integrate(p_scad, a = 3.7, lambda = 0, 0, i)[[1]]/i
}

lines(1:100, integral_vals, col = "green")

