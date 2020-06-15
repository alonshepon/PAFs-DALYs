### DOSE-RESPONSE META-ANALYSIS

## required packages
library(bd)
library(dosresmeta)
library(rms)
library(readxl)
library(dbplyr)
library(MASS)
library(ggplot2)



## Relative risk
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

## plot dose-response model
plot_dosres <-
  function(fit, df) {
    xmax <- 10 * ceiling(max(df$dose) / 10)
    
    newdata = data.frame(dose = seq(0, xmax, 10))
    with(predict(fit, newdata, xref = 0, exp = TRUE),{
      plot(get("rcs(dose, knots)dose"),
           pred, type = "l", ylab = "Relative risk", las = 1,
           xlab = "dose", ylim = c(0.5, 1.5), bty = "l")
      lines(get("rcs(dose, knots)dose"), ci.lb, lty = "dashed")
      lines(get("rcs(dose, knots)dose"), ci.ub, lty = "dashed")
    })
    rug(df$dose, quiet = T)
    points(df$dose, df$rr,
           cex = scale(df$peryears, center = FALSE), col = "blue")
  }

plot_dosres(fit_crc,crc)

#------consumption
dta <- read.csv("consumption.csv")
hist(dta$x, breaks = 20)


#compute PAF - bootstrap method
rr_fun <- rcsplineFunction(attr(fit_crc$model[[2]], "parms"), coef(fit_crc))
PRR <- mean(exp(rr_fun(dta$x)))    #exponent is to convert from log
#in this computation P is assumed to be equal to 1/n n number of points
(PRR - 1) / PRR


#Method of moments > PAF
#coefficents of Gamma distribution
m <- mean(dta$x)
v <- var(dta$x)
b <- m / v
a <- m * b
#
qqplot(dgamma((dta$x), a, b),(dta$x))
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
(PRR1 - 1) / PRR1


#fitting Gamma by maximum likelihood

dta2 <- dta$x
dta2[dta2 == 0] <- 1e-2
fit_mle <- fitdistr(dta2, "gamma")
gam=data.frame(y=dgamma(dta$x, fit_mle$estimate[1], fit_mle$estimate[2]),x=dta$x);

qqplot(dgamma(dta$x, fit_mle$estimate[1], fit_mle$estimate[2]) ,(dta$x))
ggplot(data=dta, aes(x=x)) + 
  #geom_histogram(aes(x=x)) +
geom_point(data=gam,aes(y=y,x=x))
  int <-integrate(
    function(x)
      dgamma(x, fit_mle$estimate[1], fit_mle$estimate[2]) *
      exp(rr_fun(x)),
    lower = 0,
    upper = Inf)
PRR2 <- int$value
(PRR2 - 1) / PRR2



