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

####################
## Galaxy dataset ##
####################

## Load galaxy dataset.

data("galaxy")

galaxyest <- list(normal = NULL, lognormal = NULL, Weibull = NULL, gamma = NULL)

## Estimate number of components, component weights and component parameters.

pdf <- c("normal", "lognormal", "Weibull", "gamma")

for (i in 1:4) {
  galaxyest[[i]] <- REBMIX(Dataset = list(galaxy = galaxy),
    Preprocessing = c("histogram", "Parzen window"),
    cmax = 10,
    Criterion = c("AIC", "BIC"),
    Variables = "continuous",
    pdf = pdf[i],
    K = 7:20)
}

summary(galaxyest$normal)
summary(galaxyest$lognormal)
summary(galaxyest$Weibull)
summary(galaxyest$gamma)

plot(galaxyest$Weibull, pos = 1, what = c("den", "dis"), ncol = 2, npts = 1000)

coef(galaxyest$Weibull, pos = 1)
