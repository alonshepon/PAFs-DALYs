### DOSE-RESPONSE META-ANALYSIS

## required packages
library(bd)
library(dosresmeta)
library(rms)

## import data
crc <- readxl("non-linear-crc.xlsx")

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
