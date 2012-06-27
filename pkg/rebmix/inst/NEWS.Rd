\name{NEWS} \title{News for Package \pkg{rebmix}}
\newcommand{\cpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

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
\section{Changes in version 2.4.1}{ 
\itemize{ 
\item Continuous parametric family \code{"gamma"} is added.  
}
}


