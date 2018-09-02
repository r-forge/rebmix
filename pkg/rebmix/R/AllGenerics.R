setGeneric("RNGMIX", 
  function(model = "RNGMIX",
    Dataset.name = character(),
    rseed = -1,
    n = numeric(),
    Theta = list(), ...)
  standardGeneric("RNGMIX"))

setGeneric("REBMIX", 
  function(model = "REBMIX",
    Dataset = list(),
    Preprocessing = character(),
    cmax = 15,
    Criterion = "AIC",
    pdf = character(),
    theta1 = numeric(),
    theta2 = numeric(),
    K = "auto",
    y0 = numeric(),
    ymin = numeric(),
    ymax = numeric(),
    ar = 0.1,
    Restraints = "loose", ...)
  standardGeneric("REBMIX"))
  
setGeneric("coef", 
  function(x = NULL,
    pos = 1, ...)
  standardGeneric("coef"))
  
setGeneric(".IC", function(x = NULL, Criterion = "AIC", pos = 1, ...) standardGeneric(".IC")) 
  
setGeneric("logL", function(x = NULL, pos = 1, ...) standardGeneric("logL")) 
setGeneric("AIC", function(x = NULL, pos = 1, ...) standardGeneric("AIC"))     
setGeneric("AIC3", function(x = NULL, pos = 1, ...) standardGeneric("AIC3"))    
setGeneric("AIC4", function(x = NULL, pos = 1, ...) standardGeneric("AIC4"))    
setGeneric("AICc", function(x = NULL, pos = 1, ...) standardGeneric("AICc"))    
setGeneric("BIC", function(x = NULL, pos = 1, ...) standardGeneric("BIC"))    
setGeneric("CAIC", function(x = NULL, pos = 1, ...) standardGeneric("CAIC"))    
setGeneric("HQC", function(x = NULL, pos = 1, ...) standardGeneric("HQC"))    
setGeneric("MDL2", function(x = NULL, pos = 1, ...) standardGeneric("MDL2"))    
setGeneric("MDL5", function(x = NULL, pos = 1, ...) standardGeneric("MDL5"))    
setGeneric("AWE", function(x = NULL, pos = 1, ...) standardGeneric("AWE"))    
setGeneric("CLC", function(x = NULL, pos = 1, ...) standardGeneric("CLC"))   
setGeneric("ICL", function(x = NULL, pos = 1, ...) standardGeneric("ICL"))   
setGeneric("ICLBIC", function(x = NULL, pos = 1, ...) standardGeneric("ICLBIC"))   
setGeneric("PRD", function(x = NULL, pos = 1, ...) standardGeneric("PRD"))   
setGeneric("SSE", function(x = NULL, pos = 1, ...) standardGeneric("SSE"))   
setGeneric("PC", function(x = NULL, pos = 1, ...) standardGeneric("PC"))

setGeneric("Zp", function(x = NULL, s = expression(c)) standardGeneric("Zp"))

setGeneric("demix", 
  function(x = NULL, 
    pos = 1, 
    variables = expression(1:d), ...)
  standardGeneric("demix"))

setGeneric("pemix", 
  function(x = NULL,
    pos = 1, 
    variables = expression(1:d),
    lower.tail = TRUE, 
    log.p = FALSE, ...) 
  standardGeneric("pemix"))
  
setGeneric("dfmix",  
  function(x = NULL,
    Dataset = data.frame(),
    pos = 1, 
    variables = expression(1:d), ...)
  standardGeneric("dfmix"))
  
setGeneric("pfmix",  
  function(x = NULL,
    Dataset = data.frame(),  
    pos = 1, 
    variables = expression(1:d),
    lower.tail = TRUE, 
    log.p = FALSE, ...) 
  standardGeneric("pfmix"))
  
setGeneric("split", 
  function(p = NULL,
    Dataset = data.frame(),  
    class = numeric(), ...) 
  standardGeneric("split"))
  
setGeneric("chunk", 
  function(x = NULL,
    variables = expression(1:d)) 
  standardGeneric("chunk"))    

setGeneric("boot", 
  function(x = NULL,
    rseed = -1,  
    pos = 1,
    Bootstrap = "parametric",
    B = 100, 
    n = numeric(), 
    replace = TRUE,
    prob = numeric(), ...)
  standardGeneric("boot"))
  
setGeneric("RCLRMIX", 
  function(model = "RCLRMIX",
    x = NULL,
    pos = 1, 
    Zt = factor(), ...)
  standardGeneric("RCLRMIX"))

setGeneric("RCLSMIX", 
  function(model = "RCLSMIX",
    x = list(),
    Dataset = data.frame(), 
    Zt = factor(), ...)
  standardGeneric("RCLSMIX"))
  
setGeneric("BFSMIX",
  function(model = "RCLSMIX",
    x = list(),
    Dataset = data.frame(),
    Zt = factor(), ...)
  standardGeneric("BFSMIX"))
  
setGeneric("g.d", function(x = NULL) standardGeneric("g.d"))
setGeneric("g.Theta1<-", function(x, i, value) standardGeneric("g.Theta1<-"))
setGeneric("g.Theta2<-", function(x, i, value) standardGeneric("g.Theta2<-"))
  
setGeneric("g.Dataset.name", function(x = NULL) standardGeneric("g.Dataset.name"))
setGeneric("g.rseed", function(x = NULL) standardGeneric("g.rseed"))
setGeneric("g.n", function(x = NULL) standardGeneric("g.n"))
setGeneric("g.Theta", function(x = NULL, pos = 0) standardGeneric("g.Theta"))
setGeneric("g.Dataset", function(x = NULL, pos = 0, ...) standardGeneric("g.Dataset"))
setGeneric("g.Zt", function(x = NULL) standardGeneric("g.Zt"))
setGeneric("g.w", function(x = NULL, pos = 0) standardGeneric("g.w"))
setGeneric("g.Variables", function(x = NULL) standardGeneric("g.Variables"))
setGeneric("g.ymin", function(x = NULL) standardGeneric("g.ymin"))
setGeneric("g.ymax", function(x = NULL) standardGeneric("g.ymax"))

setGeneric("g.Preprocessing", function(x = NULL) standardGeneric("g.Preprocessing"))
setGeneric("g.cmax", function(x = NULL) standardGeneric("g.cmax"))
setGeneric("g.Criterion", function(x = NULL) standardGeneric("g.Criterion"))
setGeneric("g.pdf", function(x = NULL) standardGeneric("g.pdf"))
setGeneric("g.theta1", function(x = NULL) standardGeneric("g.theta1"))
setGeneric("g.theta2", function(x = NULL) standardGeneric("g.theta2"))
setGeneric("g.K", function(x = NULL) standardGeneric("g.K"))
setGeneric("g.y0", function(x = NULL) standardGeneric("g.y0"))
setGeneric("g.ar", function(x = NULL) standardGeneric("g.ar"))
setGeneric("g.Restraints", function(x = NULL) standardGeneric("g.Restraints"))
setGeneric("g.summary", function(x = NULL) standardGeneric("g.summary"))
setGeneric("g.pos", function(x = NULL) standardGeneric("g.pos"))
setGeneric("g.opt.c", function(x = NULL) standardGeneric("g.opt.c"))
setGeneric("g.opt.IC", function(x = NULL) standardGeneric("g.opt.IC"))
setGeneric("g.opt.logL", function(x = NULL) standardGeneric("g.opt.logL"))
setGeneric("g.opt.D", function(x = NULL) standardGeneric("g.opt.D"))
setGeneric("g.all.K", function(x = NULL) standardGeneric("g.all.K"))
setGeneric("g.all.IC", function(x = NULL) standardGeneric("g.all.IC"))

setGeneric("g.Bootstrap", function(x = NULL) standardGeneric("g.Bootstrap"))
setGeneric("g.B", function(x = NULL) standardGeneric("g.B"))
setGeneric("g.replace", function(x = NULL) standardGeneric("g.replace"))
setGeneric("g.prob", function(x = NULL) standardGeneric("g.prob"))
setGeneric("g.c", function(x = NULL) standardGeneric("g.c"))
setGeneric("g.c.se", function(x = NULL) standardGeneric("g.c.se"))
setGeneric("g.c.cv", function(x = NULL) standardGeneric("g.c.cv"))
setGeneric("g.c.mode", function(x = NULL) standardGeneric("g.c.mode"))
setGeneric("g.c.prob", function(x = NULL) standardGeneric("g.c.prob"))
setGeneric("g.w.se", function(x = NULL) standardGeneric("g.w.se"))
setGeneric("g.w.cv", function(x = NULL) standardGeneric("g.w.cv"))
setGeneric("g.Theta.se", function(x = NULL) standardGeneric("g.Theta.se"))
setGeneric("g.Theta.cv", function(x = NULL) standardGeneric("g.Theta.cv"))

setGeneric("g.Zp", function(x = NULL) standardGeneric("g.Zp"))
setGeneric("g.p", function(x = NULL) standardGeneric("g.p"))
setGeneric("g.pi", function(x = NULL) standardGeneric("g.pi"))
setGeneric("g.P", function(x = NULL) standardGeneric("g.P"))
setGeneric("g.tau", function(x = NULL) standardGeneric("g.tau"))
setGeneric("g.from", function(x = NULL) standardGeneric("g.from"))
setGeneric("g.to", function(x = NULL) standardGeneric("g.to"))
setGeneric("g.EN", function(x = NULL) standardGeneric("g.EN"))
setGeneric("g.ED", function(x = NULL) standardGeneric("g.ED"))

setGeneric("g.o", function(x = NULL) standardGeneric("g.o"))
setGeneric("g.s", function(x = NULL) standardGeneric("g.s"))
setGeneric("g.ntrain", function(x = NULL) standardGeneric("g.ntrain"))
setGeneric("g.ntest", function(x = NULL) standardGeneric("g.ntest"))
setGeneric("g.CM", function(x = NULL) standardGeneric("g.CM"))
setGeneric("g.Accuracy", function(x = NULL) standardGeneric("g.Accuracy"))
setGeneric("g.Error", function(x = NULL) standardGeneric("g.Error"))
setGeneric("g.Precision", function(x = NULL) standardGeneric("g.Precision"))
setGeneric("g.Sensitivity", function(x = NULL) standardGeneric("g.Sensitivity"))
setGeneric("g.Specificity", function(x = NULL) standardGeneric("g.Specificity"))
setGeneric("g.Chunks", function(x = NULL) standardGeneric("g.Chunks"))

setGeneric("g.levels", function(x = NULL) standardGeneric("g.levels"))
setGeneric("g.train", function(x = NULL) standardGeneric("g.train"))
setGeneric("g.test", function(x = NULL) standardGeneric("g.test"))







