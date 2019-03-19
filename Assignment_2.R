### Dr. Cao - Assignment 2 ###
### Barinder Thind ###
### Functional Data Analysis ###

### Loading Libraries
library(tidyverse)
library(fda)
library(refund)
library(fda.usc)
library(ggfortify)
library(cluster)
###

# provide working directory of the data file
Gene <- read.table("https://raw.githubusercontent.com/caojiguo/FDAworkshop/master/Assignments/gene58rep34.txt", header = T)

# Looking at data
head(Gene)
unique(Gene$gene)
str(Gene)
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

## First, let's create the functional observations

# creating bspline basis
bsplinebasis_a2 <- create.bspline.basis(rangeval = c(1, 10),
                                        norder = 4, 
                                        nbasis = 5)

# Making functional observations
gene_fda = smooth.basis(timepts,
                        as.matrix(t(gene.10[,2:11])), 
                        bsplinebasis_a2)

# Plotting 
plot(gene_fda$fd, 
     xlab='Time Point', 
     ylab='Gene Expression')

##### (a) - Doing the fPCA

# Running fPCA
gene_pca = pca.fd(gene_fda$fd, nharm = 4)

# Plotting
plot(gene_pca$values)

gene_pca$scores

# Getting fPCAs
gene_fd_pca = gene_pca$harmonics
gene_fd_pca_vals = eval.fd(timepts, gene_fd_pca)

# Plotting all the fPCAs
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
matplot(timepts, gene_fd_pca_vals, xlab = 'Time Point', ylab='PCs',
        lwd=4,lty=1, cex.lab=2.5, cex.axis=2.5, type='l')
legend(8.2,-0.2,c('PC1','PC2','PC3','PC4'),col=1:4,lty=1,lwd=5)
title('Gene Principle Component Functions')

### (b)

# Looking at proportion of variance explained by the PC's
gene_pca$varprop

# Clustering first 2 PCs since they explain nearly 92% of the variance and
# also, it is easier to visualise the data this way
pca_scores <- data.frame(gene_pca$scores, 
                         cluster = as.factor(pam(gene_pca$scores[,1:2], 4)$clustering))

# Plotting the PC's with clusters
pca_scores %>% 
  ggplot(aes(x = X1, y = X2, color = cluster)) + 
  geom_point() +
  stat_ellipse(aes(x=X1, y=X2,color=cluster), type = "norm", level = 0.85) +
  theme_bw() +
  labs(title = "Clustering First 2 PC's", x = "PC 1", y = "PC 2") + 
  theme(plot.title = element_text(hjust = 0.5))

# Clustered plot of original data
data.frame(gene.10, cluster = pca_scores$cluster) %>%
  mutate(obs_in_exp = row_number()) %>%
  gather(time, gene_expression, contains("time")) %>%
  mutate(time_point = parse_number(time),
         cluster = as.factor(cluster)) %>%
  ggplot(aes(x = time_point, 
             y = gene_expression, 
             group = obs_in_exp, 
             colour = cluster)) +
  geom_line(size = 0.5) +
  labs(x = "Time", y = "Gene Expression Value", title = "Clustered Data") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# Easier to look at using facet-wrap
data.frame(gene.10, cluster = pca_scores$cluster) %>%
  mutate(obs_in_exp = row_number()) %>%
  gather(time, gene_expression, contains("time")) %>%
  mutate(time_point = parse_number(time),
         cluster = as.factor(cluster)) %>%
  ggplot(aes(x = time_point, 
             y = gene_expression, 
             group = obs_in_exp, 
             colour = cluster)) +
  geom_line(size = 0.5) +
  facet_wrap(~cluster) +
  labs(x = "Time", y = "Gene Expression Value", title = "Clustered Data") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

