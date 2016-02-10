setMethod("plot", 
          signature(x = "RCLSMVNORM", y = "missing"),
function(x,
  y,
  pos = 1,
  what = c("density"),
  nrow = 1,
  ncol = 1,
  npts = 200,
  n = 200,
  cex = 0.8,
  fg = "black",
  lty = "solid",
  lwd = 1,
  pty = "m",
  tcl = 0.5,
  plot.cex = 0.8,
  plot.pch = 19,
  contour.drawlabels = FALSE,
  contour.labcex = 0.8,
  contour.method = "flattest",
  contour.nlevels = 12, ...)
{
  if (missing(x)) {
    stop(sQuote("x"), " object of class RCLSMVNORM is requested!", call. = FALSE)
  }

  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(pos) <- 1

  if ((pos < 1) || (pos > nrow(x@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(x@summary), "!", call. = FALSE)
  }
  
  if (!is.character(what)) {
    stop(sQuote("what"), " character vector is requested!", call. = FALSE)
  }

  what <- match.arg(what, .rebmix.plot$what, several.ok = TRUE)  

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

  if (!is.wholenumber(npts)) {
    stop(sQuote("npts"), " integer is requested!", call. = FALSE)
  }

  if (npts < 1) {
    stop(sQuote("npts"), " must be greater than 0!", call. = FALSE)
  }
  
  if (!is.wholenumber(n)) {
    stop(sQuote("n"), " integer is requested!", call. = FALSE)
  }

  if (n < 1) {
    stop(sQuote("n"), " must be greater than 0!", call. = FALSE)
  }

  d <- nrow(x@Dataset)

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

#  if (.Device == "tikz output") {
#    item <- list()
#
#    item[[1]] <- "$\\textrm{Dataset}$"
#    item[[2]] <- "$\\; = \\;$"
#    item[[3]] <- paste("$\\textrm{", x@summary[pos, "Dataset"], "}$, ", sep = "")
#
#    item[[4]] <- "$\\textrm{Preprocessing}$"
#    item[[5]] <- "$\\; = \\;$"
#    item[[6]] <- paste("$\\textrm{", x@summary[pos, "Preprocessing"], "}$, ", sep = "")
#
#    item[[7]] <- "$\\textrm{Restraints}$"
#    item[[8]] <- "$\\; = \\;$"
#    item[[9]] <- paste("$\\textrm{", x@summary[pos, "Restraints"], "}$, ", sep = "")
#
#    item[[10]] <- "$D$"
#    item[[11]] <- ""
#    item[[12]] <- ""
#
#    item[[13]] <- "$c_{\\mathrm{max}}$"
#    item[[14]] <- "$\\; = \\;$"
#    item[[15]] <- paste("$", x@summary[pos, "cmax"], "$, ", sep = "")
#
#    item[[16]] <- "$a_{\\mathrm{r}}$"
#    item[[17]] <- "$\\; = \\;$"
#    item[[18]] <- paste("$", as.number(x@summary[pos, "ar"]), "$, ", sep = "")
#
#    item[[19]] <- "$c$"
#    item[[20]] <- "$\\; = \\;$"
#    item[[21]] <- paste("$", x@summary[pos, "c"], "$, ", sep = "")
#
#    item[[22]] <- ""
#    item[[23]] <- ""
#    item[[24]] <- ""
#
#    item[[28]] <- paste("$\\mathrm{", x@summary[pos, "Criterion"], "}$", sep = "")
#    item[[29]] <- "$\\; = \\;$"
#    item[[30]] <- paste("$", as.number(x@summary[pos, "IC"]), "$, ", sep = "")
#
#    item[[31]] <- "$\\mathrm{log}\\, L$"
#    item[[32]] <- "$\\; = \\;$"
#    item[[33]] <- paste("$", as.number(x@summary[pos, "logL"]), "$.", sep = "")
#
#    i <- 1; legend <- list(); legend[[i]] <- item[[1]]
#
#    for (j in c(2:9, 13:21, 25:33)) {
#      legendwidth <- strwidth(paste(legend[[i]], item[[j]], sep = ""), units = "figure", cex = 1.0)
#
#      if (legendwidth > ncol) {
#        i <- i + 1; legend[[i]] <- item[[j]]
#      }
#      else {
#        legend[[i]] <- paste(legend[[i]], item[[j]], sep = "")
#      }
#    }
#  }
#  else {
#    item <- list()
#
#    item[[1]] <- "Dataset"
#    item[[2]] <- " = "
#    item[[3]] <- paste(x@summary[pos, "Dataset"], ", ", sep = "")
#
#    item[[4]] <- "Preprocessing"
#    item[[5]] <- " = "
#    item[[6]] <- paste(x@summary[pos, "Preprocessing"], ", ", sep = "")
#
#    item[[7]] <- "Restraints"
#    item[[8]] <- " = "
#    item[[9]] <- paste(x@summary[pos, "Restraints"], ", ", sep = "")
#
#    item[[10]] <- "D"
#    item[[11]] <- ""
#    item[[12]] <- ""
#
#    item[[13]] <- bquote(c[max])
#    item[[14]] <- " = "
#    item[[15]] <- paste(x@summary[pos, "cmax"], ", ", sep = "")
#
#    item[[16]] <- bquote(a[r])
#    item[[17]] <- " = "
#    item[[18]] <- paste(as.number(x@summary[pos, "ar"]), ", ", sep = "")
#
#    item[[19]] <- "c"
#    item[[20]] <- " = "
#    item[[21]] <- paste(x@summary[pos, "c"], ", ", sep = "")
#
#    item[[22]] <- ""
#    item[[23]] <- ""
#    item[[24]] <- ""
#
#    item[[28]] <- as.character(x@summary[pos, "Criterion"])
#    item[[29]] <- " = "
#    item[[30]] <- paste(as.number(x@summary[pos, "IC"]), ", ", sep = "")
#
#    item[[31]] <- "log L"
#    item[[32]] <- " = "
#    item[[33]] <- paste(as.number(x@summary[pos, "logL"]), ".", sep = "")
#
#    i <- 1; legend <- list(); legend[[i]] <- bquote(.(item[[1]]))
#
#    for (j in c(2:9, 13:21, 25:33)) {
#      legendwidth <- strwidth(bquote(paste(.(legend[[i]]), .(item[[j]]), sep = "")), units = "figure", cex = 1.0)
#
#      if (legendwidth > ncol) {
#        i <- i + 1; legend[[i]] <- item[[j]]
#      }
#      else {
#        legend[[i]] <- bquote(paste(.(legend[[i]]), .(item[[j]]), sep = ""))
#      }
#    }
#  }

  legend <- list("ABCD")
  
  par(oma = c(length(legend) + 0.2, 0.2, 0.2, 0.2))

  ey <- as.matrix(x@Dataset)

  if (N > 0) {
      figno <- 0
    
      ramp <- colorRamp(colors = c("magenta", "blue", "cyan", "green", "yellow", "red"),
        space = "rgb",
        interpolate = "linear")

      for (i in 1:(d - 1)) {
        for (j in (i + 1):d) {
          plot(x = ey[, i],
            y = ey[, j],
            type = "p",
            main = "",
            sub = "",
            xlab = "",
            ylab = "",
            col = clRed,
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
            for (l in 1:length(legend)) {
              mtext(text = legend[[l]],
                side = 1,
                line = l - 1,
                outer = TRUE,
                adj = 0.5,
                padj = 0.2,
                cex = cex)
            }

            figno <- 0
          }
        }
      }
    }
  }
  
  rm(list = ls()[!(ls() %in% c("opar"))])

  invisible(opar)
}) # plot         
