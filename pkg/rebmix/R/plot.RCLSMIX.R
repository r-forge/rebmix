setMethod("plot", 
          signature(x = "RCLSMIX", y = "missing"),
function(x,
  y,
  nrow = 1,
  ncol = 1,
  cex = 0.8,
  fg = "black",
  lty = "solid",
  lwd = 1,
  pty = "m",
  tcl = 0.5,
  plot.cex = 0.8,
  plot.pch = 19, ...)
{
  if (missing(x)) {
    stop(sQuote("x"), " object of class RCLSMIX is requested!", call. = FALSE)
  }

  if (!is.wholenumber(nrow)) {
    stop(sQuote("nrow"), " integer is requested!", call. = FALSE)
  }

  if (nrow < 1) {
    stop(sQuote("nrow"), " must be greater than 0!", call. = FALSE)
  }

  if (!is.wholenumber(ncol)) {
    stop(sQuote("ncol"), " integer is requested!", call. = FALSE)
  }

  if (ncol < 1) {
    stop(sQuote("ncol"), " must be greater than 0!", call. = FALSE)
  }

  d <- ncol(x@Dataset)

  nrow <- max(1, nrow)
  ncol <- max(1, ncol)

  N <- d * (d - 1) / 2
  
  opar <- par(mfrow = c(nrow, ncol),
    cex = cex,
    cex.axis = 1.0,
    fg = fg,
    lty = lty,
    lwd = lwd,
    mar = c(1.2, 1.2, 1.2, 1.2),
    oma = c(1.2, 0.2, 0.2, 0.2),
    pty = pty,
    tcl = tcl, ...)
  
  par(oma = c(length("X") + 0.2, 0.2, 0.2, 0.2))

  ey <- as.matrix(x@Dataset); et <- as.numeric(x@Zt) - 1; ep <- as.numeric(x@Zp) - 1
  
  if (N > 0) {
    figno <- 0
    
    ramp <- colorRamp(colors = c("magenta", "blue", "cyan", "green", "yellow", "red"),
      space = "rgb",
      interpolate = "linear")
        
    zlim <- c(0, x@s - 1); zmax <- zlim[2]
      
    plot.col <- rgb(ramp(ep / zmax), maxColorValue = 255); plot.col[et != ep] <- "black" 
    
    legend.col <- c(rgb(ramp(zlim[1]:zlim[2] / zmax), maxColorValue = 255), "black")

    if (.Device == "tikz output") {
      legend <- c(paste("$", 1:x@s, "$", sep = ""), "$\\mathrm{Error}$")
    }
    else {
      legend <- c(paste(bquote(.(1:x@s)), sep = ""), "Error")  
    }   

    for (i in 1:(d - 1)) {
      for (j in (i + 1):d) {
        plot(x = ey[, i],
          y = ey[, j],
          type = "p",
          main = "",
          sub = "",
          xlab = "",
          ylab = "",
          col = plot.col,
          axes = FALSE,
          lwd = 1,
          cex = plot.cex,
          pch = plot.pch)

        box(col = fg, lty = "solid", lwd = 1)

        axis(side = 3,
          outer = FALSE,
          lty = "solid",
          lwd = 1,
          hadj = 0.5,
          padj = 1.0)

        axis(side = 2,
          outer = FALSE,
          lty = "solid",
          lwd = 1,
          hadj = 0.5,
          padj = 1.0)

        if (.Device == "tikz output") {
          text <- paste("$y_{", i, "}$", "$\\; - \\;$", "$y_{", j, "}$", sep = "")
        }
        else {
          text <- bquote(y[.(i)] - y[.(j)])
        }

        mtext(text = text,
          side = 1,
          line = 0,
          outer = FALSE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)
          
        figno <- figno + 1

        if ((figno == nrow * ncol) || ((i == d - 1) && (j == d))) {
          par(fig = c(0, 1, 0, 1), 
            oma = c(0, 0, 0, 0), 
            mar = c(0, 0, 0, 0), 
            new = TRUE)
          
          plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
          
          legend("bottom", 
            legend = legend,
            col = legend.col,
            lty = 0,
            pch = plot.pch,
            bty = "n",
            cex = 1.0,
            horiz = TRUE,
            inset = c(0, 0),
            xpd = TRUE)
  
          par(mfrow = c(nrow, ncol),
            cex = cex,
            cex.axis = 1.0,
            fg = fg,
            lty = lty,
            lwd = lwd,
            mar = c(1.2, 1.2, 1.2, 1.2),
            oma = c(1.2, 0.2, 0.2, 0.2),
            pty = pty,
            tcl = tcl, ...)
  
          par(oma = c(length("X") + 0.2, 0.2, 0.2, 0.2))

          figno <- 0
        }
      }
    }
  }
  
  rm(list = ls()[!(ls() %in% c("opar"))])

  invisible(opar)
}) # plot         

setMethod("plot", 
          signature(x = "RCLSMVNORM", y = "missing"),
function(x,
  y,
  nrow = 1,
  ncol = 1,
  cex = 0.8,
  fg = "black",
  lty = "solid",
  lwd = 1,
  pty = "m",
  tcl = 0.5,
  plot.cex = 0.8,
  plot.pch = 19, ...)
{
  if (missing(x)) {
    stop(sQuote("x"), " object of class RCLSMVNORM is requested!", call. = FALSE)
  }

  if (!is.wholenumber(nrow)) {
    stop(sQuote("nrow"), " integer is requested!", call. = FALSE)
  }

  if (nrow < 1) {
    stop(sQuote("nrow"), " must be greater than 0!", call. = FALSE)
  }

  if (!is.wholenumber(ncol)) {
    stop(sQuote("ncol"), " integer is requested!", call. = FALSE)
  }

  if (ncol < 1) {
    stop(sQuote("ncol"), " must be greater than 0!", call. = FALSE)
  }

  d <- ncol(x@Dataset)

  nrow <- max(1, nrow)
  ncol <- max(1, ncol)

  N <- d * (d - 1) / 2
  
  opar <- par(mfrow = c(nrow, ncol),
    cex = cex,
    cex.axis = 1.0,
    fg = fg,
    lty = lty,
    lwd = lwd,
    mar = c(1.2, 1.2, 1.2, 1.2),
    oma = c(1.2, 0.2, 0.2, 0.2),
    pty = pty,
    tcl = tcl, ...)
  
  par(oma = c(length("X") + 0.2, 0.2, 0.2, 0.2))

  ey <- as.matrix(x@Dataset); et <- as.numeric(x@Zt) - 1; ep <- as.numeric(x@Zp) - 1
  
  if (N > 0) {
    figno <- 0
    
    ramp <- colorRamp(colors = c("magenta", "blue", "cyan", "green", "yellow", "red"),
      space = "rgb",
      interpolate = "linear")
        
    zlim <- c(0, x@s - 1); zmax <- zlim[2]
      
    plot.col <- rgb(ramp(ep / zmax), maxColorValue = 255); plot.col[et != ep] <- "black" 
    
    legend.col <- c(rgb(ramp(zlim[1]:zlim[2] / zmax), maxColorValue = 255), "black")
          
    if (.Device == "tikz output") {
      legend <- c(paste("$", 1:x@s, "$", sep = ""), "$\\mathrm{Error}$")
    }
    else {
      legend <- c(paste(bquote(.(1:x@s)), sep = ""), "Error")  
    }

    for (i in 1:(d - 1)) {
      for (j in (i + 1):d) {
        plot(x = ey[, i],
          y = ey[, j],
          type = "p",
          main = "",
          sub = "",
          xlab = "",
          ylab = "",
          col = plot.col,
          axes = FALSE,
          lwd = 1,
          cex = plot.cex,
          pch = plot.pch)

        box(col = fg, lty = "solid", lwd = 1)

        axis(side = 3,
          outer = FALSE,
          lty = "solid",
          lwd = 1,
          hadj = 0.5,
          padj = 1.0)

        axis(side = 2,
          outer = FALSE,
          lty = "solid",
          lwd = 1,
          hadj = 0.5,
          padj = 1.0)

        if (.Device == "tikz output") {
          text <- paste("$y_{", i, "}$", "$\\; - \\;$", "$y_{", j, "}$", sep = "")
        }
        else {
          text <- bquote(y[.(i)] - y[.(j)])
        }

        mtext(text = text,
          side = 1,
          line = 0,
          outer = FALSE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)
          
        figno <- figno + 1

        if ((figno == nrow * ncol) || ((i == d - 1) && (j == d))) {
          par(fig = c(0, 1, 0, 1), 
            oma = c(0, 0, 0, 0), 
            mar = c(0, 0, 0, 0), 
            new = TRUE)
          
          plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
          
          legend("bottom", 
            legend = legend,
            col = legend.col,
            lty = 0,
            pch = plot.pch,
            bty = "n",
            cex = 1.0,
            horiz = TRUE,
            inset = c(0, 0),
            xpd = TRUE)
  
          par(mfrow = c(nrow, ncol),
            cex = cex,
            cex.axis = 1.0,
            fg = fg,
            lty = lty,
            lwd = lwd,
            mar = c(1.2, 1.2, 1.2, 1.2),
            oma = c(1.2, 0.2, 0.2, 0.2),
            pty = pty,
            tcl = tcl, ...)
  
          par(oma = c(length("X") + 0.2, 0.2, 0.2, 0.2))

          figno <- 0
        }
      }
    }
  }
  
  rm(list = ls()[!(ls() %in% c("opar"))])

  invisible(opar)
}) # plot         
