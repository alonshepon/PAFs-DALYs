library("data.table")
library("reshape2")
library("ggplot2")
library("dplyr")
library("tidyverse")


## import data
crc <- read_xlsx("non-linear-crc.xlsx")

## fit dose-response model
fit_dosres <-
  function(x) {
    knots <- quantile(x$dose, c(.1, .5, .9))
    fit <-
      dosresmeta(formula = logrr ~ rcs(dose, knots),
                 type = "ir",
                 id = factor(study),
                 se = se,
                 cases = cases,
                 n = peryears,
                 data = x,
                 method = "fixed")
  }
fit_crc <- fit_dosres(crc)
rr_fun <- rcsplineFunction(attr(fit_crc$model[[2]], "parms"), coef(fit_crc))


#------load consumption data
dta <- read.csv("consumption.csv")
leg=100;
R <- data.frame(matrix(ncol = 3, nrow = 100))
x <- c("mean", "std", "PAF")
colnames(R) <- x
d<-1
for (val in seq(from = 0.5,to = 2,length.out = 10)){
  for (val1 in seq(from = 0.5,to = 2,length.out = 10)){
    #browser()
    #Method of moments > PAF
    #coefficents of Gamma distribution
    m <- mean(dta$x)*val
    v <- var(dta$x)*val1
    b <- m / v
    a <- m * b
    #If X is a gamma(α, β) random variable and the shape parameter α is large relative to the scale parameter β, 
    #then X approximately has a normal random variable with the same mean and variance.
    
    
    y<-rgamma(leg, a, b)  #randomly generate gamma distribution
    op=fitdist(y, "gamma", method = "mme")  #fit curve
    int <-integrate(
      function(x)
        dgamma(x, coef(op)[1], coef(op)[2]) *
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


## single graph
#coefficents of Gamma distribution

fg_try <- fitdist(dta$x, "gamma",method = "mme")
plot(fg_try)
    