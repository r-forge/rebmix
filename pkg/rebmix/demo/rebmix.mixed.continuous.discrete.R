##############################################
## R sources for reproducing the results in ##
##   Marko Nagode:                          ##
##   rebmix: An R Package for Continuous    ##
##   and Discrete Finite Mixture Models     ##
##############################################

options(prompt = "> ", continue = "+ ", width = 70,
  useFancyQuotes = FALSE, digits = 3)

###################
## Preliminaries ##
###################

## load package.

library("rebmix")

#######################################
## Mixed continuous-discrete dataset ##
#######################################

## Generate mixed continuous-discrete dataset.

n <- c(400, 100, 500)

Theta <- rbind(pdf1 = rep("lognormal", 3),
  theta1.1 = c(1.0, 3.5, 2.5), theta2.1 = c(0.3, 0.2, 0.4),
  pdf2 = rep("Poisson", 3),
  theta1.2 = c(2.0, 10.0, 15.0),
  pdf3 = rep("binomial", 3),
  theta1.3 = c(10, 10, 10), theta2.3 = c(0.9, 0.1, 0.7),
  pdf4 = rep("Weibull", 3),
  theta1.4 = c(2.0, 10.0, 25.0), theta2.4 = c(3.0, 7.0, 20.0))

mixed <- RNGMIX(Dataset = "mixed", n = n, Theta = Theta)

## Estimate number of components, component weights and component parameters.

Sturges <- as.integer(1 + log2(sum(n))) ## Minimum v follows the Sturges rule.
Log10 <- as.integer(10 * log10(sum(n))) ## Maximum v follows the Log10 rule.

mixedest <- REBMIX(Dataset = mixed$Dataset, Preprocessing = "histogram", cmax = 9,
  Criterion = "BIC", Variables = c("continuous", "discrete", "discrete", "continuous"),
  pdf = c("lognormal", "Poisson", "binomial", "Weibull"),
  Theta1 = c(NA, NA, 10, NA), K = kseq(Sturges, Log10, 0.1))

plot(mixedest, what = c("dens", "marg", "IC", "D"), nrow = 4, ncol = 3, npts = 1000)

## Bootstrap finite mixture.

mixedboot <- boot.REBMIX(x = mixedest, pos = 1, Bootstrap = "p", 
  B = 100, n = NULL, replace = TRUE, prob = NULL)

mixedboot
