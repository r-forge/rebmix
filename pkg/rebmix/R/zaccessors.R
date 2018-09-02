setMethod("g.c", signature(x = "THETA"), function(x) x@c)
setMethod("g.d", signature(x = "THETA"), function(x) x@d)
setMethod("g.pdf", signature(x = "THETA"), function(x) x@pdf)
setMethod("g.Theta", signature(x = "THETA"), function(x) x@Theta)

setMethod("g.Theta1<-", 
  signature(x = "THETA"),
function(x, i, value)
{
  x@Theta[[2 + (i - 1) * 3]] <- value
})

setMethod("g.Theta2<-", 
  signature(x = "THETA"),
function(x, i, value)
{
  x@Theta[[3 + (i - 1) * 3]] <- value
})

setMethod("g.Dataset.name", signature(x = "RNGMIX"), function(x) x@Dataset.name)
setMethod("g.rseed", signature(x = "RNGMIX"), function(x) x@rseed)
setMethod("g.n", signature(x = "RNGMIX"), function(x) x@n)
setMethod("g.Theta", signature(x = "RNGMIX"), function(x) x@Theta)

setMethod("g.Dataset", 
          signature(x = "RNGMIX"), 
function(x, pos)
{
  if ((pos < 1) || (pos > length(x@Dataset))) {
    output <- x@Dataset
  }
  else {
    output <- x@Dataset[[pos]]
  }
  
  output
}) ## Dataset
          
setMethod("g.Zt", signature(x = "RNGMIX"), function(x) x@Zt)
setMethod("g.w", signature(x = "RNGMIX"), function(x) x@w)
setMethod("g.Variables", signature(x = "RNGMIX"), function(x) x@Variables)
setMethod("g.ymin", signature(x = "RNGMIX"), function(x) x@ymin)
setMethod("g.ymax", signature(x = "RNGMIX"), function(x) x@ymax)

setMethod("g.Dataset", 
          signature(x = "REBMIX"), 
function(x, pos)
{
  if ((pos < 1) || (pos > length(x@Dataset))) {
    output <- x@Dataset
  }
  else {
    output <- x@Dataset[[pos]]
  }
  
  output
}) ## Dataset

setMethod("g.Preprocessing", signature(x = "REBMIX"), function(x) x@Preprocessing)
setMethod("g.cmax", signature(x = "REBMIX"), function(x) x@cmax)
setMethod("g.Criterion", signature(x = "REBMIX"), function(x) x@Criterion)
setMethod("g.Variables", signature(x = "REBMIX"), function(x) x@Variables)
setMethod("g.pdf", signature(x = "REBMIX"), function(x) x@pdf)
setMethod("g.theta1", signature(x = "REBMIX"), function(x) x@theta1)
setMethod("g.theta2", signature(x = "REBMIX"), function(x) x@theta2)
setMethod("g.K", signature(x = "REBMIX"), function(x) x@K)
setMethod("g.y0", signature(x = "REBMIX"), function(x) x@y0)
setMethod("g.ymin", signature(x = "REBMIX"), function(x) x@ymin)
setMethod("g.ymax", signature(x = "REBMIX"), function(x) x@ymax)
setMethod("g.ar", signature(x = "REBMIX"), function(x) x@ar)
setMethod("g.Restraints", signature(x = "REBMIX"), function(x) x@Restraints)

setMethod("g.w", 
          signature(x = "REBMIX"), 
function(x, pos)
{
  if ((pos < 1) || (pos > length(x@w))) {
    output <- x@w
  }
  else {
    output <- x@w[[pos]]
  }
  
  output
}) ## w

setMethod("g.Theta", 
          signature(x = "REBMIX"), 
function(x, pos)
{
  if ((pos < 1) || (pos > length(x@Theta))) {
    output <- x@Theta
  }
  else {
    output <- x@Theta[[pos]]
  }
  
  output
}) ## Theta

setMethod("g.summary", signature(x = "REBMIX"), function(x) x@summary)
setMethod("g.pos", signature(x = "REBMIX"), function(x) x@pos)
setMethod("g.opt.c", signature(x = "REBMIX"), function(x) x@opt.c)
setMethod("g.opt.IC", signature(x = "REBMIX"), function(x) x@opt.IC)
setMethod("g.opt.logL", signature(x = "REBMIX"), function(x) x@opt.logL)
setMethod("g.opt.D", signature(x = "REBMIX"), function(x) x@opt.D)
setMethod("g.all.K", signature(x = "REBMIX"), function(x) x@all.K)
setMethod("g.all.IC", signature(x = "REBMIX"), function(x) x@all.IC)

setMethod("g.rseed", signature(x = "REBMIX.boot"), function(x) x@rseed)
setMethod("g.pos", signature(x = "REBMIX.boot"), function(x) x@pos)
setMethod("g.Bootstrap", signature(x = "REBMIX.boot"), function(x) x@Bootstrap)
setMethod("g.B", signature(x = "REBMIX.boot"), function(x) x@B)
setMethod("g.n", signature(x = "REBMIX.boot"), function(x) x@n)
setMethod("g.replace", signature(x = "REBMIX.boot"), function(x) x@replace)
setMethod("g.prob", signature(x = "REBMIX.boot"), function(x) x@prob)
setMethod("g.c", signature(x = "REBMIX.boot"), function(x) x@c)
setMethod("g.c.se", signature(x = "REBMIX.boot"), function(x) x@c.se)
setMethod("g.c.cv", signature(x = "REBMIX.boot"), function(x) x@c.cv)
setMethod("g.c.mode", signature(x = "REBMIX.boot"), function(x) x@c.mode)
setMethod("g.c.prob", signature(x = "REBMIX.boot"), function(x) x@c.prob)
setMethod("g.w", signature(x = "REBMIX.boot"), function(x) x@w)
setMethod("g.w.se", signature(x = "REBMIX.boot"), function(x) x@w.se)
setMethod("g.w.cv", signature(x = "REBMIX.boot"), function(x) x@w.cv)
setMethod("g.Theta", signature(x = "REBMIX.boot"), function(x) x@Theta)
setMethod("g.Theta.se", signature(x = "REBMIX.boot"), function(x) x@Theta.se)
setMethod("g.Theta.cv", signature(x = "REBMIX.boot"), function(x) x@Theta.cv)

setMethod("g.pos", signature(x = "RCLRMIX"), function(x) x@pos)
setMethod("g.Zt", signature(x = "RCLRMIX"), function(x) x@Zt)
setMethod("g.Zp", signature(x = "RCLRMIX"), function(x) x@Zp)
setMethod("g.c", signature(x = "RCLRMIX"), function(x) x@c)
setMethod("g.p", signature(x = "RCLRMIX"), function(x) x@p)
setMethod("g.pi", signature(x = "RCLRMIX"), function(x) x@pi)
setMethod("g.P", signature(x = "RCLRMIX"), function(x) x@P)
setMethod("g.tau", signature(x = "RCLRMIX"), function(x) x@tau)
setMethod("g.prob", signature(x = "RCLRMIX"), function(x) x@prob)
setMethod("g.from", signature(x = "RCLRMIX"), function(x) x@from)
setMethod("g.to", signature(x = "RCLRMIX"), function(x) x@to)
setMethod("g.EN", signature(x = "RCLRMIX"), function(x) x@EN)
setMethod("g.ED", signature(x = "RCLRMIX"), function(x) x@ED)

setMethod("g.s", signature(x = "RCLS.chunk"), function(x) x@s)
setMethod("g.levels", signature(x = "RCLS.chunk"), function(x) x@levels)
setMethod("g.ntrain", signature(x = "RCLS.chunk"), function(x) x@ntrain)
setMethod("g.train", signature(x = "RCLS.chunk"), function(x) x@train)
setMethod("g.ntest", signature(x = "RCLS.chunk"), function(x) x@ntest)
setMethod("g.test", signature(x = "RCLS.chunk"), function(x) x@test)
setMethod("g.Zt", signature(x = "RCLS.chunk"), function(x) x@Zt)

setMethod("g.o", signature(x = "RCLSMIX"), function(x) x@o)
setMethod("g.Dataset", signature(x = "RCLSMIX"), function(x) x@Dataset)
setMethod("g.s", signature(x = "RCLSMIX"), function(x) x@s)
setMethod("g.ntrain", signature(x = "RCLSMIX"), function(x) x@ntrain)
setMethod("g.P", signature(x = "RCLSMIX"), function(x) x@P)
setMethod("g.ntest", signature(x = "RCLSMIX"), function(x) x@ntest)
setMethod("g.Zt", signature(x = "RCLSMIX"), function(x) x@Zt)
setMethod("g.Zp", signature(x = "RCLSMIX"), function(x) x@Zp)
setMethod("g.CM", signature(x = "RCLSMIX"), function(x) x@CM)
setMethod("g.Accuracy", signature(x = "RCLSMIX"), function(x) x@Accuracy)
setMethod("g.Error", signature(x = "RCLSMIX"), function(x) x@Error)
setMethod("g.Precision", signature(x = "RCLSMIX"), function(x) x@Precision)
setMethod("g.Sensitivity", signature(x = "RCLSMIX"), function(x) x@Sensitivity)
setMethod("g.Specificity", signature(x = "RCLSMIX"), function(x) x@Specificity)
setMethod("g.Chunks", signature(x = "RCLSMIX"), function(x) x@Chunks)
