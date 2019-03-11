### Dr. Cao - Assignment 2 ###
### Barinder Thind ###
### Functional Data Analysis ###

### Loading Libraries
library(tidyverse)
library(fda)
library(refund)
library(fda.usc)
###

# provide working directory of the data file
Gene <- read.table("https://raw.githubusercontent.com/caojiguo/FDAworkshop/master/Assignments/gene58rep34.txt", header = T)

# Looking at data
head(Gene)
unique(Gene$gene)
table(Gene$gene)
# =============================================

# take F10 out
gene.10 <- Gene[Gene$gene == "F10",]
head(gene.10)
dim(gene.10)

# observation time points
timepts = seq(1 ,10 ,1)
# open a new graph window in R in a Windows laptop .
# run quartz () if you are using a Mac computer
windows()
matplot(t(gene.10[,2:11]), type ="l", xlab = "Time point", ylab = "Gene expression")


#### (Preprocessing)

## First, let's create the functional observations using a roughness penalty
# creating bspline basis
bsplinebasis_a2 <- create.bspline.basis(rangeval = c(1, 10),
                                        norder = 4, 
                                        nbasis = 5)

# Choosing Smoothing Parameters using cross-validation
lambda_val = seq(0, 1, 0.01)
len_lam = length(lambda_val)
gcv_val = c()

# Running loop
for (i in 1:len_lam) {
  
  # Particular value of lambda
  lambda = lambda_val[i]
  
  # Setting up
  fdParobj = fdPar(bsplinebasis_a2, 2, lambda)
  
  # Creating FD Object
  smoothlist = smooth.basis(timepts, as.matrix(t(gene.10[,2:11])),
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

# Making functional observations
# Creating best FD Object
gene_fda_rough = smooth.basisPar(timepts,
                                    as.matrix(t(gene.10[,2:11])), 
                                    bsplinebasis_a2, 
                                    lambda = lambda_val[which.min(gcv_val)])

# Plotting 
plot(tempfd$fd,xlab='day',ylab='temperature',cex.lab=1.5,cex.axis=1.5)

##### (a) - Doing the fPCA

# Running fPCA
gene_pca = pca.fd(gene_fda_rough$fd, nharm = 4)

# Plotting
plot(gene_pca$values)

# Getting fPCAs
gene_fd_pca = gene_pca$harmonics
geen_fd_pca_vals = eval.fd(timepts, gene_fd_pca)

# Plotting all the fPCAs
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
matplot(timepts, geen_fd_pca_vals, xlab='Time Point',ylab='PCs',
        lwd=4,lty=1,cex.lab=2.5,cex.axis=2.5,type='l')
legend(0,-0.07,c('PC1','PC2','PC3','PC4'),col=1:4,lty=1,lwd=5)
title('Gene Principle Component Functions')

### (b)


