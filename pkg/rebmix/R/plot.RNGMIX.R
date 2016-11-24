setMethod("plot", 
          signature(x = "RNGMIX", y = "missing"),
function(x,
  y,
  pos = 1,  
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
    stop(sQuote("x"), " object of class RNGMIX is requested!", call. = FALSE)
  }
  
  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(pos) <- 1

  if ((pos < 1) || (pos > length(x@Dataset))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", length(x@Dataset), "!", call. = FALSE)
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

  Zt <- as.numeric(levels(x@Zt))[x@Zt]
  
  zlim <- c(0, max(1, max(Zt) - 1)); zmax <- zlim[2] 

  d <- ncol(x@Dataset[[pos]])

  nrow <- max(1, nrow)
  ncol <- max(1, ncol)

  N <- d * (d - 1) / 2
  
  opar <- list(); ipar <- 1  
  
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

  par(oma = c(1 + 0.2, 0.2, 0.2, 0.2))
  
  ey <- as.matrix(x@Dataset[[pos]]); et <- Zt - 1
  
  ramp <- colorRamp(colors = c("magenta", "blue", "cyan", "green", "yellow", "red"),
    space = "rgb",
    interpolate = "linear")
    
  unique.Zt <- unique(Zt); s <- length(unique.Zt); sort.unique.Zt <- sort(unique.Zt)
  
  plot.col <- rgb(ramp(et / zmax), maxColorValue = 255)
  
  legend.col <- rgb(ramp((sort.unique.Zt - 1) / zmax), maxColorValue = 255)
  
  legend.text <- as.character(sort.unique.Zt)
  
  legend.pch <- rep(plot.pch, s)  
  
  if (N > 0) {
    figno <- 0

    for (i in 1:(d - 1)) {
      for (j in (i + 1):d) {
        plot(x = ey[, i],
          y = ey[, j],
          type = "p",
          main = "",
          sub = "",
          xlab = "",
          ylab = "",
          xlim = range(ey[, i]),
          ylim = range(ey[, j]),
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
          
          .legendA(s = s, text = legend.text, col = legend.col, pch = legend.pch, error = 0)          
  
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
    
          par(oma = c(1 + 0.2, 0.2, 0.2, 0.2))

          figno <- 0
        }
        
        opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1        
      }
    }
  }
  else {
    plot(x = ey[, 1],
      y = et + 1,
      type = "p",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      xlim = range(ey[, 1]),
      ylim = range(et + 1),      
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
      text <- paste("$y_{1}$", "$\\; - \\;$", "$Z_{t}(y_{1})$", sep = "")
    }
    else {
      text <- bquote(y[1] - Z[t](y[1]))
    }

    mtext(text = text,
      side = 1,
      line = 0,
      outer = FALSE,
      adj = 0.5,
      padj = 0.2,
      cex = cex)

    par(fig = c(0, 1, 0, 1), 
      oma = c(0, 0, 0, 0), 
      mar = c(0, 0, 0, 0), 
      new = TRUE)
          
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    
    .legendA(s = s, text = legend.text, col = legend.col, pch = legend.pch, error = 0)          
      
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
      
    par(oma = c(1 + 0.2, 0.2, 0.2, 0.2))
    
    opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1    
  }
  
  rm(list = ls()[!(ls() %in% c("opar"))])

  invisible(opar)
}) # plot

setMethod("plot", 
          signature(x = "RNGMVNORM", y = "missing"),
function(x,
  y,
  pos = 1,  
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
    stop(sQuote("x"), " object of class RNGMVNORM is requested!", call. = FALSE)
  }

  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(pos) <- 1

  if ((pos < 1) || (pos > length(x@Dataset))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", length(x@Dataset), "!", call. = FALSE)
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

  Zt <- as.numeric(levels(x@Zt))[x@Zt]
  
  zlim <- c(0, max(1, max(Zt) - 1)); zmax <- zlim[2] 

  d <- ncol(x@Dataset[[pos]])

  nrow <- max(1, nrow)
  ncol <- max(1, ncol)

  N <- d * (d - 1) / 2
  
  opar <- list(); ipar <- 1  
  
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

  par(oma = c(1 + 0.2, 0.2, 0.2, 0.2))
  
  ey <- as.matrix(x@Dataset[[pos]]); et <- Zt - 1
  
  ramp <- colorRamp(colors = c("magenta", "blue", "cyan", "green", "yellow", "red"),
    space = "rgb",
    interpolate = "linear")
    
  unique.Zt <- unique(Zt); s <- length(unique.Zt); sort.unique.Zt <- sort(unique.Zt)
  
  plot.col <- rgb(ramp(et / zmax), maxColorValue = 255)
  
  legend.col <- rgb(ramp((sort.unique.Zt - 1) / zmax), maxColorValue = 255)
  
  legend.text <- as.character(sort.unique.Zt)
  
  legend.pch <- rep(plot.pch, s)  
  
  if (N > 0) {
    figno <- 0

    for (i in 1:(d - 1)) {
      for (j in (i + 1):d) {
        plot(x = ey[, i],
          y = ey[, j],
          type = "p",
          main = "",
          sub = "",
          xlab = "",
          ylab = "",
          xlim = range(ey[, i]),
          ylim = range(ey[, j]),
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
          
          .legendA(s = s, text = legend.text, col = legend.col, pch = legend.pch, error = 0)          
  
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
    
          par(oma = c(1 + 0.2, 0.2, 0.2, 0.2))

          figno <- 0
        }
        
        opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1        
      }
    }
  }
  else {
    plot(x = ey[, 1],
      y = et + 1,
      type = "p",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      xlim = range(ey[, 1]),
      ylim = range(et + 1),      
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
      text <- paste("$y_{1}$", "$\\; - \\;$", "$Z_{t}(y_{1})$", sep = "")
    }
    else {
      text <- bquote(y[1] - Z[t](y[1]))
    }

    mtext(text = text,
      side = 1,
      line = 0,
      outer = FALSE,
      adj = 0.5,
      padj = 0.2,
      cex = cex)

    par(fig = c(0, 1, 0, 1), 
      oma = c(0, 0, 0, 0), 
      mar = c(0, 0, 0, 0), 
      new = TRUE)
          
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    
    .legendA(s = s, text = legend.text, col = legend.col, pch = legend.pch, error = 0)          
      
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
      
    par(oma = c(1 + 0.2, 0.2, 0.2, 0.2))
    
    opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1    
  }
  
  rm(list = ls()[!(ls() %in% c("opar"))])

  invisible(opar)
}) # plot
