\name{NEWS} \title{News for Package \pkg{rebmix}}
\newcommand{\cpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}
\section{Changes in version 2.6.1}{ 
\itemize{
\item Renewed vignette. 
\item Renewed references in Reference manual.
\item Method \code{PD} is added.  
\item Fixed bug in C functions \code{REBMIXKNN()}, \code{REBMIXPW()} and \code{REBMIXH()}.
}
}
\section{Changes in version 2.6.0}{ 
\itemize{ 
\item Switch \code{_TIME_LEFT_SWITCH} introduced in C code.
\item Simplified treatment of loose restraints in C functions \code{RoughEstimationKNN()}, \code{RoughEstimationPW()} and \code{RoughEstimationH()}.
\item Fixed bug in C functions \code{REBMIXKNN()}, \code{REBMIXPW()} and \code{REBMIXH()}.
\item Minimum weight multiplier \code{b} is eliminated.
\item Fixed source and destination overlap in memcpy bug in C code.
}
}
\section{Changes in version 2.5.1}{ 
\itemize{ 
\item Fixed checking top-level files.
\item Fixed checking for old-style vignette sources.
\item Packages mixtools and mclust excluded from vignette.
}
}
\section{Changes in version 2.5.0}{ 
\itemize{ 
\item Rough parameter estimation is improved by introducing constrained entropy maximization.
\item Fixed bug in C functions \code{RoughEstimationKNN()}, \code{RoughEstimationPW()} and \code{RoughEstimationH()}.
\item Fixed bug in \code{boot.REBMIX} method.
}
}
\section{Changes in version 2.4.2}{ 
\itemize{ 
\item Fixed symbol eating bug for Mac OS X.  
}
}
\section{Changes in version 2.4.1}{ 
\itemize{ 
\item Continuous parametric family \code{"gamma"} is added.  
}
}
\section{Changes in version 2.4.0}{ 
\itemize{ 
\item Methods \code{AIC}, \code{AIC3}, \code{AIC4}, \code{AICc} and \code{CAIC} are added. 
\item Method \code{AWE} is added. 
\item Method \code{BIC} is added. 
\item Method \code{boot.REBMIX} for bootstrap is added. 
\item Method \code{CLC} is added. 
\item Function \code{demix} for empirical density calculation is added. 
\item Function \code{dfmix} for predictive density calculation is added. 
\item Method \code{HQC} is added. 
\item Method \code{ICL} is added. 
\item Method \code{ICLBIC} is added. 
\item Method \code{logL} is added. 
\item Methods \code{MDL2} and \code{MDL5} are added. 
\item Function \code{pemix} for empirical distribution function calculation is added. 
\item Function \code{pfmix} for predictive distribution function calculation is added. 
\item Function \code{plot.REBMIX} is extended. It enables besides \code{"density"} plots also \code{"marginal"}, 
\code{"c"}, \code{"IC"}, \code{"logL"}, \code{"D"} and \code{"distribution"} plots.
\item Additional arguments \code{what} enters method \code{plot.REBMIX}.
\item Method \code{PRD} is added. 
\item Methods \code{print.boot.REBMIX} and \code{summary.boot.REBMIX} are added.
\item Method \code{SSE} is added. 
\item Dataset \code{weibull} is added. 
}
}
\section{Changes in version 2.3.0}{ 
\itemize{ 
\item Discrete parametric families \code{"binomial"}, \code{"Poisson"} and \code{"Dirac"} are added. 
\item Information criteria \code{"D"} and \code{"SSE"} are added.
\item Methods \code{print.RNGMIX}, \code{print.REBMIX}, \code{coef.REBMIX} and \code{summary.REBMIX} are added.
\item Method \code{predict.list} for class membership prediction is added.
\item The package does not write to stdout and read from stdout any more.
\item The package is protected against wrong input arguments.
\item Additional arguments \code{Variables}, \code{Theta1}, \code{Theta2}, \code{ymin}, \code{ymax} and \code{b} 
    enter function \code{REBMIX}.
\item Minimum radius \code{RMIN} is eliminated.
\item Simplified treatment of global modes with equal probabilities in C functions \code{GlobalModeKNN()}, 
    \code{GlobalModePW()} and \code{GlobalModeH()}.
\item Improved treatment of potentially singular solutions in C functions \code{EnhancedEstimationKNN()}, 
    \code{EnhancedEstimationPW()} and \code{EnhancedEstimationH()}.
\item Improved convergence criteria in C functions \code{REBMIXKNN()}, \code{REBMIXPW()} and \code{REBMIXH()}.
\item Fixed bug in C function \code{PreprocessingPW()}.
\item Fixed bug in C functions \code{RoughEstimationKNN()}, \code{RoughEstimationPW()} and 
    \code{RoughEstimationH()}.
\item Fixed bug in \code{plot.REBMIX} method. 
\item Dataset \code{adult} is added. 
}
}
