setGeneric("RNGMIX",
  function(model = NULL,
    Dataset.name = character(0),
    rseed = -1,
    n = numeric(0),
    Theta = list(0))
  standardGeneric("RNGMIX"))

setGeneric("REBMIX",
  function(model = NULL,
    Dataset = list(0),
    Preprocessing = character(0),
    cmax = 15,
    Criterion = "AIC",
    pdf = character(0),
    theta1 = numeric(0),
    theta2 = numeric(0),
    K = numeric(0),
    y0 = numeric(0),
    ymin = numeric(0),
    ymax = numeric(0),
    ar = 0.1,
    Restraints = "loose", ...)
  standardGeneric("REBMIX"))

