\name{mmenusoptd.mae}
\alias{mmenusoptd.mae}
\title{
Creates the main menu tcltk window
}
\description{
Alternative to directly using the function \code{\link{soptdmaeA}} on the R console, the function \code{mmenusoptd.mae} creates the main menu tcltk window with the list of four optimality criteria (A, MV, D and E) for search of sequential optimal or near-optimal block or row-column designs for two-colour cDNA microarray experiments. }
\usage{
mmenusoptd.mae()
}
\details{
The main menu window created using the function \code{mmenusoptd.mae} contains the list of four optimality criteria (the A-,
MV -, D- and E) for sequential block or row-column designs that allows user to choose an optimality criterion of interest. The function 
will then call for a function \code{fixparsoptd.mae} for further option regarding the setup of parametric combinations, initial optimal or near-optimal design, and output, see \code{fixparsoptd.mae} documentation. 
}
\author{
Dibaba Bayisa Gemechu, Legesse Kassa Debusho, and Linda Haines
}

\seealso{
\code{\link{soptdmaeA}}, \code{\link{fixparsoptd.mae}}, \code{\link{tcltk}}, \code{\link{TkWidgets}}
}
