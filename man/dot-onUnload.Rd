% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\name{.onUnload}
\alias{.onUnload}
\title{Package Unload Hook}
\usage{
.onUnload(libpath)
}
\arguments{
\item{libpath}{Character string giving the library where the package is installed}
}
\description{
Performs cleanup operations when the package is unloaded.
Ensures proper cleanup of parallel processing and package environment.
}
\keyword{internal}
