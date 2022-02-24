setMethod("plot",
          signature(x = "REBMIX", y = "missing"),
function(x,
  y,
  pos = 1,
  what = c("pdf"),
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
  contour.nlevels = 12, 
  log = "", ...)
{
  if (missing(x)) {
    stop(sQuote("x"), " object of class REBMIX is requested!", call. = FALSE)
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

  ni <- ncol(x@summary)

  Theta <- .extractThetaA(x@w[[pos]], x@Theta[[pos]])

  d <- nrow(Theta)

  nrow <- max(1, nrow)
  ncol <- max(1, ncol)

  N <- 0

  opar <- list(); ipar <- 1
  
  opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1

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

  C <- x@summary[pos, "Preprocessing"]

  if (is.na(C)) {
    C <- "NA"
    
    items <- c(2:3, 19:21, 31:33)
  }
  else {
    items <- c(2:6, 13:21, 25:33)
  }

  item <- list()

  item[[1]] <- "Dataset"
  item[[2]] <- " = "
  item[[3]] <- paste(x@summary[pos, "Dataset"], ", ", sep = "")

  item[[4]] <- "Preprocessing"
  item[[5]] <- " = "
  item[[6]] <- paste(x@summary[pos, "Preprocessing"], ", ", sep = "")

  item[[7]] <- ""
  item[[8]] <- ""
  item[[9]] <- ""

  item[[10]] <- ""
  item[[11]] <- ""
  item[[12]] <- ""

  item[[13]] <- bquote(c[max])
  item[[14]] <- " = "
  item[[15]] <- paste(x@summary[pos, "cmax"], ", ", sep = "")

  item[[16]] <- bquote(a[r])
  item[[17]] <- " = "
  item[[18]] <- paste(as.number(x@summary[pos, "ar"]), ", ", sep = "")

  item[[19]] <- "c"
  item[[20]] <- " = "
  item[[21]] <- paste(x@summary[pos, "c"], ", ", sep = "")

  item[[22]] <- ""
  item[[23]] <- ""
  item[[24]] <- ""

  if (C == .rebmix$Preprocessing[1]) {
    item[[25]] <- "v"
    item[[26]] <- " = "
    item[[27]] <- paste(x@summary[pos, "v/k"], ", ", sep = "")
  }
  else
  if (C == .rebmix$Preprocessing[2]) {
    item[[25]] <- "v"
    item[[26]] <- " = "
    item[[27]] <- paste(x@summary[pos, "v/k"], ", ", sep = "")
  }
  else
  if (C == .rebmix$Preprocessing[3]) {
    item[[25]] <- "k"
    item[[26]] <- " = "
    item[[27]] <- paste(x@summary[pos, "v/k"], ", ", sep = "")
  }

  item[[28]] <- as.character(x@summary[pos, "Criterion"])
  item[[29]] <- " = "
  item[[30]] <- paste(as.number(x@summary[pos, "IC"]), ", ", sep = "")

  item[[31]] <- "log L"
  item[[32]] <- " = "
  item[[33]] <- paste(as.number(x@summary[pos, "logL"]), ".", sep = "")

  i <- 1; legend <- list(); legend[[i]] <- bquote(.(item[[1]]))

  for (j in items) {
    legendwidth <- strwidth(bquote(paste(.(legend[[i]]), .(item[[j]]), sep = "")), units = "figure", cex = 1.0)

    if (legendwidth > ncol) {
      i <- i + 1; legend[[i]] <- item[[j]]
    }
    else {
      legend[[i]] <- bquote(paste(.(legend[[i]]), .(item[[j]]), sep = ""))
    }
  }

  par(oma = c(length(legend) + 0.2, 0.2, 0.2, 0.2))

# Dataset <- as.character(x@summary[pos, "Dataset"])
  
  Dataset <- x@Dataset[[which(names(x@Dataset) == x@summary[pos, "Dataset"])]]
  
  if (class(Dataset) == "data.frame") {
    Y.type <- 0
    
    ey <- as.matrix(Dataset)
  }
  else
  if (class(Dataset) == "Histogram") {
    Y.type <- 1
    
    ey <- as.matrix(Dataset@Y)
  }
  
# ey <- as.matrix(x@Dataset[[which(names(x@Dataset) == x@summary[pos, "Dataset"])]])

  y0 <- array(data = 0.0, dim = d, dimnames = NULL)
  h <- array(data = 0.0, dim = d, dimnames = NULL)

  lim <- array(data = 0.0, dim = c(2, d), dimnames = NULL)

  py <- list(d)

  Variables <- match.arg(x@Variables, .rebmix$Variables, several.ok = TRUE)

  pdf <- match.arg(x@pdf, .rebmix$pdf, several.ok = TRUE)

  for (i in 1:d) {
    if (Y.type == 0) {
      if (C == .rebmix$Preprocessing[1]) {
        k <- as.numeric(x@summary[pos, "v/k"])
        y0[i] <- as.numeric(x@summary[pos, paste("y0", if (d > 1) i, sep = "")])
        h[i] <- as.numeric(x@summary[pos, paste("h", if (d > 1) i, sep = "")])

        lim[, i] <- range(ey[, i], finite = TRUE)
      }
      else
      if (C == .rebmix$Preprocessing[2]) {
        h[i] <- as.numeric(x@summary[pos, paste("h", if (d > 1) i, sep = "")])

        lim[, i] <- range(ey[, i], finite = TRUE)
      }
      else
      if (C == .rebmix$Preprocessing[3]) {
        k <- as.numeric(x@summary[pos, "v/k"])

        h[i] <- as.numeric(x@summary[pos, paste("h", if (d > 1) i, sep = "")])

        lim[, i] <- range(ey[, i], finite = TRUE)
      }
      else {
        lim[, i] <- range(ey[, i], finite = TRUE)
      }
    }
    else
    if (Y.type == 1) {
      v <- nrow(ey)
      
      h[i] <- Dataset@h[i]
    
      lim[, i] <- range(ey[, i], finite = TRUE)
    
    }

    if (abs(lim[2, i] - lim[1, i]) < 1e-6) {
      lim[2, i] <- lim[1, i] + 1.0
    }

    if (Variables[i] == .rebmix$Variables[2]) {
      py[[i]] <- sort(unique(ey[, i]))
    }
    else {
      py[[i]] <- seq(from = lim[1, i], to = lim[2, i], length.out = npts)
    }
  }

  w <- as.numeric(x@w[[pos]])

  if (d > 1) {
    if (any(match(.rebmix.plot$what[1], what, nomatch = 0))) {
      N <- d * (d - 1) / 2

      figno <- 0

      ramp <- colorRamp(colors = c("magenta", "blue", "cyan", "green", "yellow", "red"),
        space = "rgb",
        interpolate = "linear")

      for (i in 1:(d - 1)) {
        for (j in (i + 1):d) {
          pdens <- outer(py[[i]], py[[j]], ".dfmix.xy", w, Theta[i,], Theta[j,])
          
          zlim <- range(pdens, finite = TRUE)
          
          if (Y.type == 0) {
            if (C == .rebmix$Preprocessing[1]) {
              edens <- .densHistogram.xy(k, ey[, i], ey[, j], y0[i], lim[, i][1], lim[, i][2], y0[j], lim[, j][1], lim[, j][2], h[i], h[j], Variables[i], Variables[j], pdf[i], pdf[j])
            }
            else
            if (C == .rebmix$Preprocessing[2]) {
              edens <- .densKDE.xy(ey[, i], ey[, j], h[i], h[j], n)
            }
            else
            if (C == .rebmix$Preprocessing[3]) {
              edens <- .densKNearestNeighbour.xy(ey[, i], ey[, j], k, h[i], h[j], n)
            }
            else {
              edens <- .densSample.xy(ey[, i], ey[, j], zlim[1], n)
            }
          }
          else
          if (Y.type == 1) {
            edens <- .densK.xy(v, ey[, i], ey[, j], ey[, d + 1], h[i], h[j])
          }

          zlim <- range(zlim, edens$z, finite = TRUE)

          plot(x = edens$x,
            y = edens$y,
            type = "p",
            main = "",
            sub = "",
            xlab = "",
            ylab = "",
            col = rgb(ramp((edens$z - zlim[1]) / (zlim[2] - zlim[1])), maxColorValue = 255),
            axes = FALSE,
            lwd = 1,
            cex = plot.cex,
            pch = plot.pch, 
            log = log)

          if ((Variables[i] == .rebmix$Variables[2]) && (Variables[j] == .rebmix$Variables[2])) {
            z <- as.vector(pdens); z <- z != 0.0

            points(x = rep(py[[i]], length(py[[j]]))[z],
              y = rep(py[[j]], each = length(py[[i]]))[z],
              type = "p",
              xlab = "",
              ylab = "",
              col = rgb(ramp((pdens[z] - zlim[1]) / (zlim[2] - zlim[1])), maxColorValue = 255),
              lwd = 1,
              cex = plot.cex * 0.5,
              pch = plot.pch)
          }
          else
          if ((Variables[i] == .rebmix$Variables[2]) && (Variables[j] == .rebmix$Variables[1])) {
            for (l in 1:length(py[[i]])) {
              tx <- rep(py[[i]][l], length(py[[j]]))
              ty <- py[[j]]

              s <- 1:(length(tx) - 1)

              segments(x0 = tx[s],
                y0 = ty[s],
                x1 = tx[s + 1],
                y1 = ty[s + 1],
                xlab = "",
                ylab = "",
                col = rgb(ramp(((pdens[l, s] + pdens[l, s + 1]) / 2.0 - zlim[1]) / (zlim[2] - zlim[1])), maxColorValue = 255),
                cex = plot.cex)
            }
          }
          else
          if ((Variables[i] == .rebmix$Variables[1]) && (Variables[j] == .rebmix$Variables[2])) {
            for (l in 1:length(py[[j]])) {
              tx <- py[[i]]
              ty <- rep(py[[j]][l], length(py[[i]]))

              s <- 1:(length(tx) - 1)

              segments(x0 = tx[s],
                y0 = ty[s],
                x1 = tx[s + 1],
                y1 = ty[s + 1],
                xlab = "",
                ylab = "",
                col = rgb(ramp(((pdens[s, l] + pdens[s + 1, l]) / 2.0 - zlim[1]) / (zlim[2] - zlim[1])), maxColorValue = 255),
                cex = plot.cex)
            }
          }
          else {
            levels <- seq(from = zlim[1], to = zlim[2], length.out = contour.nlevels)

            contour(x = py[[i]],
              y = py[[j]],
              z = pdens,
              levels = levels,
              xlim = lim[, i],
              ylim = lim[, j],
              zlim = zlim,
              labcex = contour.labcex, drawlabels = contour.drawlabels, method = contour.method,
              axes = FALSE, frame.plot = FALSE,
              col = rgb(ramp((levels - zlim[1]) / (zlim[2] - zlim[1])), maxColorValue = 255),
              add = TRUE)
          }

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

          text <- bquote(y[.(i)] - y[.(j)])

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

          opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1
        }
      }
    }
  }

  m <- nrow * ncol * ceiling(N / nrow / ncol) - N

  if (any(match(.rebmix.plot$what[2], what, nomatch = 0)) || ((d == 1) && any(match(.rebmix.plot$what[1], what, nomatch = 0)))) {
    for (i in 1:d) {
      pdens <- .dfmix.x(py[[i]], w, Theta[i,])
    
      ylim <- range(pdens, finite = TRUE)
    
      if (Y.type == 0) {    
        if (C == .rebmix$Preprocessing[1]) {
          edens <- .densHistogram.x(k, ey[, i], y0[i], lim[, i][1], lim[, i][2], h[i], Variables[i], pdf[i])
        }
        else
        if (C == .rebmix$Preprocessing[2]) {
          edens <- .densKDE.x(ey[, i], h[i], n)
        }
        else
        if (C == .rebmix$Preprocessing[3]) {
          edens <- .densKNearestNeighbour.x(ey[, i], k, h[i], n)
        }
        else {
          edens <- .densSample.x(ey[, i], ylim[1], n)
        }
      }
      else
      if (Y.type == 1) {
        edens <- .densK.x(v, ey[, i], ey[, d + 1], h[i])      
      }

      ylim <- range(ylim, edens$y, finite = TRUE)

      plot(x = edens$x,
        y = edens$y,
        type = "p",
        main = "",
        sub = "",
        xlab = "",
        ylab = "",
        ylim = ylim,
        col = "black",
        axes = FALSE,
        lwd = 1,
        cex = plot.cex,
        pch = plot.pch,
        log = log)

      points(x = py[[i]],
        y = pdens,
        type = "l",
        col = "black")

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

      text <- bquote(y[.(i)] - f(y[.(i)]))

      mtext(text = text,
        side = 1,
        line = 0,
        outer = FALSE,
        adj = 0.5,
        padj = 0.2,
        cex = cex)

      if (m <= 0) {
        for (l in 1:length(legend)) {
          mtext(text = legend[[l]],
            side = 1,
            line = l - 1,
            outer = TRUE,
            adj = 0.5,
            padj = 0.2,
            cex = cex)
        }

        m <- nrow * ncol - 1
      }
      else {
        m <- m - 1
      }

      opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1
    }
  }
  
  if (any(match(.rebmix.plot$what[6], what, nomatch = 0)) || ((d == 1) && any(match(.rebmix.plot$what[8], what, nomatch = 0)))) {
    for (i in 1:d) {
      pdist <- .pfmix.x(py[[i]], w, Theta[i,])
          
      if (Y.type == 0) {
        edist <- .dist.x(ey[, i], NULL, n)
      }
      else
      if (Y.type == 1) {
        edist <- .dist.x(ey[, i], ey[, d + 1], n)
      }  

      ylim <- range(pdist, edist$y, finite = TRUE)

      plot(x = edist$x,
        y = edist$y,
        type = "p",
        main = "",
        sub = "",
        xlab = "",
        ylab = "",
        ylim = ylim,
        col = "black",
        axes = FALSE,
        lwd = 1,
        cex = plot.cex,
        pch = plot.pch)

      points(x = py[[i]],
        y = pdist,
        type = "l",
        col = "black")

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

      text <- bquote(y[.(i)] - F(y[.(i)]))

      mtext(text = text,
        side = 1,
        line = 0,
        outer = FALSE,
        adj = 0.5,
        padj = 0.2,
        cex = cex)

      if (m <= 0) {
        for (l in 1:length(legend)) {
          mtext(text = legend[[l]],
            side = 1,
            line = l - 1,
            outer = TRUE,
            adj = 0.5,
            padj = 0.2,
            cex = cex)
        }

        m <- nrow * ncol - 1
      }
      else {
        m <- m - 1
      }

      opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1
    }
  }  

  if (any(match(.rebmix.plot$what[3], what, nomatch = 0)) && !is.na(x@opt.IC[[pos]][1])) {
    ylim <- range(x@opt.IC[[pos]], finite = TRUE)

    plot(x = x@opt.c[[pos]],
      y = x@opt.IC[[pos]],
      type = "o",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      ylim = ylim,
      col = "black",
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

    text <- bquote(c - .(item[[28]]))

    mtext(text = text,
      side = 1,
      line = 0,
      outer = FALSE,
      adj = 0.5,
      padj = 0.2,
      cex = cex)

    if (m <= 0) {
      for (l in 1:length(legend)) {
        mtext(text = legend[[l]],
          side = 1,
          line = l - 1,
          outer = TRUE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)
      }

      m <- nrow * ncol - 1
    }
    else {
      m <- m - 1
    }

    opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1
  }

  if (any(match(.rebmix.plot$what[4], what, nomatch = 0)) && !is.na(x@opt.logL[[pos]][1])) {
    ylim <- range(x@opt.logL[[pos]], finite = TRUE)

    plot(x = x@opt.c[[pos]],
      y = x@opt.logL[[pos]],
      type = "o",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      ylim = ylim,
      col = "black",
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

    text <- bquote(c - .(item[[31]]))

    mtext(text = text,
      side = 1,
      line = 0,
      outer = FALSE,
      adj = 0.5,
      padj = 0.2,
      cex = cex)

    if (m <= 0) {
      for (l in 1:length(legend)) {
        mtext(text = legend[[l]],
          side = 1,
          line = l - 1,
          outer = TRUE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)
      }

      m <- nrow * ncol - 1
    }
    else {
      m <- m - 1
    }

    opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1
  }

  if (any(match(.rebmix.plot$what[5], what, nomatch = 0)) && !is.na(x@opt.D[[pos]][1])) {
    ylim <- range(x@opt.D[[pos]], finite = TRUE)

    plot(x = x@opt.c[[pos]],
      y = x@opt.D[[pos]],
      type = "o",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      ylim = ylim,
      col = "black",
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

    text <- bquote(c - .(item[[10]]))

    mtext(text = text,
      side = 1,
      line = 0,
      outer = FALSE,
      adj = 0.5,
      padj = 0.2,
      cex = cex)

    if (m <= 0) {
      for (l in 1:length(legend)) {
        mtext(text = legend[[l]],
          side = 1,
          line = l - 1,
          outer = TRUE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)
      }

      m <- nrow * ncol - 1
    }
    else {
      m <- m - 1
    }

    opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1
  }

  if (any(match(.rebmix.plot$what[7], what, nomatch = 0)) && !is.na(x@all.IC[[pos]][1])) {
    ylim <- range(x@all.IC[[pos]], finite = TRUE)

    plot(x = x@all.K[[pos]],
      y = x@all.IC[[pos]],
      type = "o",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      ylim = ylim,
      col = "black",
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

    text <- bquote(K - .(item[[28]]))

    mtext(text = text,
      side = 1,
      line = 0,
      outer = FALSE,
      adj = 0.5,
      padj = 0.2,
      cex = cex)

    if (m <= 0) {
      for (l in 1:length(legend)) {
        mtext(text = legend[[l]],
          side = 1,
          line = l - 1,
          outer = TRUE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)
      }

      m <- nrow * ncol - 1
    }
    else {
      m <- m - 1
    }

    opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1
  }

  rm(list = ls()[!(ls() %in% c("opar"))])

  invisible(opar)
}) ## plot

setMethod("plot",
          signature(x = "REBMVNORM", y = "missing"),
function(x,
  y,
  pos = 1,
  what = c("pdf"),
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
  contour.nlevels = 12,
  log = "", ...)
{
  if (missing(x)) {
    stop(sQuote("x"), " object of class REBMIX is requested!", call. = FALSE)
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

  ni <- ncol(x@summary)

  Theta <- .extractThetaB(x@w[[pos]], x@Theta[[pos]])

  d <- length(Theta[[1]]$pdf)

  nrow <- max(1, nrow)
  ncol <- max(1, ncol)

  N <- 0

  opar <- list(); ipar <- 1
  
  opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1

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

  C <- x@summary[pos, "Preprocessing"]
  
  if (is.na(C)) {
    C <- "NA"
    
    items <- c(2:3, 19:21, 31:33)
  }
  else {
    items <- c(2:6, 13:21, 25:33)
  }

  item <- list()
  
  item[[1]] <- "Dataset"
  item[[2]] <- " = "
  item[[3]] <- paste(x@summary[pos, "Dataset"], ", ", sep = "")

  item[[4]] <- "Preprocessing"
  item[[5]] <- " = "
  item[[6]] <- paste(x@summary[pos, "Preprocessing"], ", ", sep = "")
  
  item[[7]] <- ""
  item[[8]] <- ""
  item[[9]] <- ""

  item[[10]] <- ""
  item[[11]] <- ""
  item[[12]] <- ""

  item[[13]] <- bquote(c[max])
  item[[14]] <- " = "
  item[[15]] <- paste(x@summary[pos, "cmax"], ", ", sep = "")

  item[[16]] <- bquote(a[r])
  item[[17]] <- " = "
  item[[18]] <- paste(as.number(x@summary[pos, "ar"]), ", ", sep = "")

  item[[19]] <- "c"
  item[[20]] <- " = "
  item[[21]] <- paste(x@summary[pos, "c"], ", ", sep = "")

  item[[22]] <- ""
  item[[23]] <- ""
  item[[24]] <- ""

  if (C == .rebmix$Preprocessing[1]) {
    item[[25]] <- "v"
    item[[26]] <- " = "
    item[[27]] <- paste(x@summary[pos, "v/k"], ", ", sep = "")
  }
  else
  if (C == .rebmix$Preprocessing[2]) {
    item[[25]] <- "v"
    item[[26]] <- " = "
    item[[27]] <- paste(x@summary[pos, "v/k"], ", ", sep = "")
  }
  else
  if (C == .rebmix$Preprocessing[3]) {
    item[[25]] <- "k"
    item[[26]] <- " = "
    item[[27]] <- paste(x@summary[pos, "v/k"], ", ", sep = "")
  }

  item[[28]] <- as.character(x@summary[pos, "Criterion"])
  item[[29]] <- " = "
  item[[30]] <- paste(as.number(x@summary[pos, "IC"]), ", ", sep = "")

  item[[31]] <- "log L"
  item[[32]] <- " = "
  item[[33]] <- paste(as.number(x@summary[pos, "logL"]), ".", sep = "")

  i <- 1; legend <- list(); legend[[i]] <- bquote(.(item[[1]]))

  for (j in items) {
    legendwidth <- strwidth(bquote(paste(.(legend[[i]]), .(item[[j]]), sep = "")), units = "figure", cex = 1.0)

    if (legendwidth > ncol) {
      i <- i + 1; legend[[i]] <- item[[j]]
    }
    else {
      legend[[i]] <- bquote(paste(.(legend[[i]]), .(item[[j]]), sep = ""))
    }
  }

  par(oma = c(length(legend) + 0.2, 0.2, 0.2, 0.2))

# Dataset <- as.character(x@summary[pos, "Dataset"])

  Dataset <- x@Dataset[[which(names(x@Dataset) == x@summary[pos, "Dataset"])]]
  
  if (class(Dataset) == "data.frame") {
    Y.type <- 0
    
    ey <- as.matrix(Dataset)
  }
  else
  if (class(Dataset) == "Histogram") {
    Y.type <- 1
    
    ey <- as.matrix(Dataset@Y)
  }

# ey <- as.matrix(x@Dataset[[which(names(x@Dataset) == x@summary[pos, "Dataset"])]])

  y0 <- array(data = 0.0, dim = d, dimnames = NULL)
  h <- array(data = 0.0, dim = d, dimnames = NULL)

  lim <- array(data = 0.0, dim = c(2, d), dimnames = NULL)

  py <- list(d)

  Variables <- match.arg(x@Variables, .rebmix$Variables, several.ok = TRUE)

  pdf <- match.arg(x@pdf, .rebmix$pdf, several.ok = TRUE)

  for (i in 1:d) {
    if (Y.type == 0) {
      if (C == .rebmix$Preprocessing[1]) {
        k <- as.numeric(x@summary[pos, "v/k"])
        y0[i] <- as.numeric(x@summary[pos, paste("y0", if (d > 1) i, sep = "")])
        h[i] <- as.numeric(x@summary[pos, paste("h", if (d > 1) i, sep = "")])

        lim[, i] <- range(ey[, i], finite = TRUE)
      }
      else
      if (C == .rebmix$Preprocessing[2]) {
        h[i] <- as.numeric(x@summary[pos, paste("h", if (d > 1) i, sep = "")])

        lim[, i] <- range(ey[, i], finite = TRUE)
      }
      else
      if (C == .rebmix$Preprocessing[3]) {
        k <- as.numeric(x@summary[pos, "v/k"])

        h[i] <- as.numeric(x@summary[pos, paste("h", if (d > 1) i, sep = "")])

        lim[, i] <- range(ey[, i], finite = TRUE)
      }
      else {
        lim[, i] <- range(ey[, i], finite = TRUE)
      }
    }
    else
    if (Y.type == 1) {
      v <- nrow(ey)
      
      h[i] <- Dataset@h[i]
    
      lim[, i] <- range(ey[, i], finite = TRUE)
    }

    if (abs(lim[2, i] - lim[1, i]) < 1e-6) {
      lim[2, i] <- lim[1, i] + 1.0
    }

    py[[i]] <- seq(from = lim[1, i], to = lim[2, i], length.out = npts)
  }

  w <- as.numeric(x@w[[pos]])

  if (d > 1) {
    if (any(match(.rebmix.plot$what[1], what, nomatch = 0))) {
      N <- d * (d - 1) / 2

      figno <- 0

      ramp <- colorRamp(colors = c("magenta", "blue", "cyan", "green", "yellow", "red"),
        space = "rgb",
        interpolate = "linear")

      for (i in 1:(d - 1)) {
        for (j in (i + 1):d) {
          pdens <- outer(py[[i]], py[[j]], ".dfmvnorm.xy", w, Theta, i, j)
          
          zlim <- range(pdens, finite = TRUE);
          
          if (Y.type == 0) {      
            if (C == .rebmix$Preprocessing[1]) {
              edens <- .densHistogram.xy(k, ey[, i], ey[, j], y0[i], lim[, i][1], lim[, i][2], y0[j], lim[, j][1], lim[, j][2], h[i], h[j], Variables[i], Variables[j], pdf[i], pdf[j])
            }
            else
            if (C == .rebmix$Preprocessing[2]) {
              edens <- .densKDE.xy(ey[, i], ey[, j], h[i], h[j], n)
            }
            else
            if (C == .rebmix$Preprocessing[3]) {
              edens <- .densKNearestNeighbour.xy(ey[, i], ey[, j], k, h[i], h[j], n)
            }
            else {
              edens <- .densSample.xy(ey[, i], ey[, j], zlim[1], n)
            }
          }
          else
          if (Y.type == 1) {
            edens <- .densK.xy(v, ey[, i], ey[, j], ey[, d + 1], h[i], h[j])
          }        
          
          zlim <- range(zlim, edens$z, finite = TRUE)

          plot(x = edens$x,
            y = edens$y,
            type = "p",
            main = "",
            sub = "",
            xlab = "",
            ylab = "",
            col = rgb(ramp((edens$z - zlim[1]) / (zlim[2] - zlim[1])), maxColorValue = 255),
            axes = FALSE,
            lwd = 1,
            cex = plot.cex,
            pch = plot.pch, 
            log = log)

          levels <- seq(from = zlim[1], to = zlim[2], length.out = contour.nlevels)
          
          contour(x = py[[i]],
            y = py[[j]],
            z = pdens,
            levels = levels,
            xlim = lim[, i],
            ylim = lim[, j],
            zlim = zlim,
            labcex = contour.labcex, drawlabels = contour.drawlabels, method = contour.method,
            axes = FALSE, frame.plot = FALSE,
            col = rgb(ramp((levels - zlim[1]) / (zlim[2] - zlim[1])), maxColorValue = 255),
            add = TRUE)

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

          text <- bquote(y[.(i)] - y[.(j)])

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

          opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1
        }
      }
    }
  }

  m <- nrow * ncol * ceiling(N / nrow / ncol) - N

  if (any(match(.rebmix.plot$what[2], what, nomatch = 0)) || ((d == 1) && any(match(.rebmix.plot$what[1], what, nomatch = 0)))) {
    for (i in 1:d) {
      pdens <- .dfmvnorm.x(py[[i]], w, Theta, i)
      
      ylim <- range(pdens, finite = TRUE)
      
      if (Y.type == 0) {      
        if (C == .rebmix$Preprocessing[1]) {
          edens <- .densHistogram.x(k, ey[, i], y0[i], lim[, i][1], lim[, i][2], h[i], Variables[i], pdf[i])
        }
        else
        if (C == .rebmix$Preprocessing[2]) {
          edens <- .densKDE.x(ey[, i], h[i], n)
        }
        else
        if (C == .rebmix$Preprocessing[3]) {
          edens <- .densKNearestNeighbour.x(ey[, i], k, h[i], n)
        }
        else {
          edens <- .densSample.x(ey[, i], ylim[1], n)
        }
      }
      else
      if (Y.type == 1) {
        edens <- .densK.x(v, ey[, i], ey[, d + 1], h[i])      
      }      
      
      ylim <- range(ylim, edens$y, finite = TRUE)

      plot(x = edens$x,
        y = edens$y,
        type = "p",
        main = "",
        sub = "",
        xlab = "",
        ylab = "",
        ylim = ylim,
        col = "black",
        axes = FALSE,
        lwd = 1,
        cex = plot.cex,
        pch = plot.pch, 
        log = log)

      points(x = py[[i]],
        y = pdens,
        type = "l",
        col = "black")

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

      text <- bquote(y[.(i)] - f(y[.(i)]))

      mtext(text = text,
        side = 1,
        line = 0,
        outer = FALSE,
        adj = 0.5,
        padj = 0.2,
        cex = cex)

      if (m <= 0) {
        for (l in 1:length(legend)) {
          mtext(text = legend[[l]],
            side = 1,
            line = l - 1,
            outer = TRUE,
            adj = 0.5,
            padj = 0.2,
            cex = cex)
        }

        m <- nrow * ncol - 1
      }
      else {
        m <- m - 1
      }

      opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1
    }
  }
  
  if (any(match(.rebmix.plot$what[6], what, nomatch = 0)) || ((d == 1) && any(match(.rebmix.plot$what[8], what, nomatch = 0)))) {
    for (i in 1:d) {
      pdist <- .pfmvnorm.x(py[[i]], w, Theta, i)
    
      if (Y.type == 0) {
        edist <- .dist.x(ey[, i], NULL, n)
      }
      else
      if (Y.type == 1) {
        edist <- .dist.x(ey[, i], ey[, d + 1], n)
      }

      ylim <- range(pdist, edist$y, finite = TRUE)    
      
      plot(x = edist$x,
        y = edist$y,
        type = "p",
        main = "",
        sub = "",
        xlab = "",
        ylab = "",
        ylim = ylim,
        col = "black",
        axes = FALSE,
        lwd = 1,
        cex = plot.cex,
        pch = plot.pch)

      points(x = py[[i]],
        y = pdist,
        type = "l",
        col = "black")

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

      text <- bquote(y[.(i)] - F(y[.(i)]))

      mtext(text = text,
        side = 1,
        line = 0,
        outer = FALSE,
        adj = 0.5,
        padj = 0.2,
        cex = cex)

      if (m <= 0) {
        for (l in 1:length(legend)) {
          mtext(text = legend[[l]],
            side = 1,
            line = l - 1,
            outer = TRUE,
            adj = 0.5,
            padj = 0.2,
            cex = cex)
        }

        m <- nrow * ncol - 1
      }
      else {
        m <- m - 1
      }

      opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1
    }
  }

  if (any(match(.rebmix.plot$what[3], what, nomatch = 0)) && !is.na(x@opt.IC[[pos]][1])) {
    ylim <- range(x@opt.IC[[pos]], finite = TRUE)

    plot(x = x@opt.c[[pos]],
      y = x@opt.IC[[pos]],
      type = "o",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      ylim = ylim,
      col = "black",
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

    text <- bquote(c - .(item[[28]]))

    mtext(text = text,
      side = 1,
      line = 0,
      outer = FALSE,
      adj = 0.5,
      padj = 0.2,
      cex = cex)

    if (m <= 0) {
      for (l in 1:length(legend)) {
        mtext(text = legend[[l]],
          side = 1,
          line = l - 1,
          outer = TRUE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)
      }

      m <- nrow * ncol - 1
    }
    else {
      m <- m - 1
    }

    opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1
  }

  if (any(match(.rebmix.plot$what[4], what, nomatch = 0)) && !is.na(x@opt.logL[[pos]][1])) {
    ylim <- range(x@opt.logL[[pos]], finite = TRUE)

    plot(x = x@opt.c[[pos]],
      y = x@opt.logL[[pos]],
      type = "o",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      ylim = ylim,
      col = "black",
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

    text <- bquote(c - .(item[[31]]))

    mtext(text = text,
      side = 1,
      line = 0,
      outer = FALSE,
      adj = 0.5,
      padj = 0.2,
      cex = cex)

    if (m <= 0) {
      for (l in 1:length(legend)) {
        mtext(text = legend[[l]],
          side = 1,
          line = l - 1,
          outer = TRUE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)
      }

      m <- nrow * ncol - 1
    }
    else {
      m <- m - 1
    }

    opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1
  }

  if (any(match(.rebmix.plot$what[5], what, nomatch = 0)) && !is.na(x@opt.D[[pos]][1])) {
    ylim <- range(x@opt.D[[pos]], finite = TRUE)

    plot(x = x@opt.c[[pos]],
      y = x@opt.D[[pos]],
      type = "o",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      ylim = ylim,
      col = "black",
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

    text <- bquote(c - .(item[[10]]))

    mtext(text = text,
      side = 1,
      line = 0,
      outer = FALSE,
      adj = 0.5,
      padj = 0.2,
      cex = cex)

    if (m <= 0) {
      for (l in 1:length(legend)) {
        mtext(text = legend[[l]],
          side = 1,
          line = l - 1,
          outer = TRUE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)
      }

      m <- nrow * ncol - 1
    }
    else {
      m <- m - 1
    }

    opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1
  }

  if (any(match(.rebmix.plot$what[7], what, nomatch = 0)) && !is.na(x@all.IC[[pos]][1])) {
    ylim <- range(x@all.IC[[pos]], finite = TRUE)

    plot(x = x@all.K[[pos]],
      y = x@all.IC[[pos]],
      type = "o",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      ylim = ylim,
      col = "black",
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

    text <- bquote(K - .(item[[28]]))

    mtext(text = text,
      side = 1,
      line = 0,
      outer = FALSE,
      adj = 0.5,
      padj = 0.2,
      cex = cex)

    if (m <= 0) {
      for (l in 1:length(legend)) {
        mtext(text = legend[[l]],
          side = 1,
          line = l - 1,
          outer = TRUE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)
      }

      m <- nrow * ncol - 1
    }
    else {
      m <- m - 1
    }

    opar[[ipar]] <- par(no.readonly = TRUE); ipar <- ipar + 1
  }

  rm(list = ls()[!(ls() %in% c("opar"))])

  invisible(opar)
}) # plot
