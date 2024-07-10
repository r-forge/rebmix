.rebmix <- structure(list(
Preprocessing = c("histogram", "kernel density estimation", "k-nearest neighbour"),
Criterion = c("AIC", "AIC3", "AIC4", "AICc", "BIC", "CAIC", "HQC", "MDL2", "MDL5", "AWE", "CLC", "ICL", "ICL-BIC", "PC", "D", "SSE"),
Variables = c("continuous", "discrete"),
pdf = c("normal", "lognormal", "Weibull", "binomial", "Poisson", "Dirac", "gamma", "uniform", "vonMises", "Gumbel"),
pdf.nargs = c(2, 2, 2, 2, 1, 1, 2, 2, 2, 3),
pdf.Variables = c("continuous", "continuous", "continuous", "discrete", "discrete", "discrete", "continuous", "continuous", "continuous", "continuous"),
Restraints = c("rigid", "loose"),
Mode = c("all", "outliers", "outliersplus"),
### Panic Branislav.
EMStrategy = c("none", "exhaustive", "best", "single"),
EMVariant = c("EM", "ECM"),
EMAcceleration = c("fixed", "line", "golden")),
### End
.Names = c("Preprocessing", "Criterion", "Variables", "pdf", "pdf.nargs", "pdf.Variables", "Restraints", "Mode", "EMStrategy", "EMVariant", "EMAcceleration"))

.rebmix.plot <- structure(list(
what = c("pdf", "marginal pdf", "IC", "logL", "D", "marginal cdf", "K", "cdf")),
.Names = c("what"))

.rebmix.boot <- structure(list(
Bootstrap = c("parametric", "nonparametric")),
.Names = c("Bootstrap"))

.optbins <- structure(list(
Rule = c("Sturges", "Log10", "RootN", "Knuth equal", "Knuth unequal")),
.Names = c("Rule"))

.rclrmix <- structure(list(
Rule = c("Entropy", "Demp")),
.Names = c("Rule"))

.error.defaults <- structure(list(
ErrorNames = c("E_OK", "E_MEM", "E_ARG", "E_CON", "E_FILE", "E_NO_SOLUTION"),
FileNames = c("base.cpp", "rngmixf.cpp", "rngmvnormf.cpp", "rebmixf.cpp", "rebmvnormf.cpp", "emf.cpp", "Rmisc.cpp", "Rrebmix.cpp", "Rrebmvnorm.cpp")),
.Names = c("ErrorNames", "FileNames"))






