##### SCAD Function #####

### Defining Function
p_scad <- function(lambda, a, beta){
  
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

  return(vals)
  
}

### Defining beta values
beta = seq(0.001, 10, 0.01)

### Runnign function
p_scad_val <- p_scad(500, 3.7, beta)

plot(beta, p_scad_val)
