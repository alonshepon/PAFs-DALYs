library("data.table")
library("reshape2")
library("ggplot2")
library("dplyr")
library("tidyverse")

R <- data.frame(matrix(ncol = 3, nrow = 25))
x <- c("mean", "std", "PAF")
colnames(df) <- x
d<-1
for (val in seq(from = 0.5,to = 3,length.out = 5)){
  for (val1 in seq(from = 0.5,to = 3,length.out = 5)){
    #browser()
    #Method of moments > PAF
    #coefficents of Gamma distribution
    m <- mean(dta$x)*val
    v <- var(dta$x)*val1
    b <- m / v
    a <- m * b
    
    args(dgamma)
    ## function (x, shape, rate = 1, scale = 1/rate, log = FALSE)
    ## NULL
    int <-integrate(
      function(x)
        dgamma(x, a, b) *
        exp(rr_fun(x)),
      lower = 0,
      upper = Inf)
    PRR1 <- int$value
    ri<-(PRR1 - 1) / PRR1
    R$PAF[d]<-ri
    R$mean[d]<-m
    R$std[d]<-v
    d<-d+1}}

##create colormap


p3<-ggplot(R, aes(x = mean, y = std, z=PAF)) +
  geom_raster(aes(fill = PAF)) +
  scale_fill_gradient(low="grey90", high="red")+
  geom_contour(aes(colour = stat(level)),show.legend = FALSE,bins=15)
p3

dta1=rgamma(100, a,b);
fg_try <- fitdist(dta1, "gamma",method = "mse")
plot(fg_try)
    