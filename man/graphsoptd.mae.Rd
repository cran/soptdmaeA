\name{graphsoptd.mae}
\alias{graphsoptd.mae}
\title{
Creates the graphical layout of resultant sequential optimal design
}
\description{
Creates the graphical layout of resultant sequential A-, MV-, D- or E-optimal or near-optimal block or row-column design on a separate pop-up GUI tcltk window.
}
\usage{
graphsoptd.mae(trt.N, blk.N, theta, soptdesF, Optcrit, strt, sary, dtype)
}
\arguments{
\item{trt.N}{
integer, specifying  number of treatments \code{v} of resultant sequential optimal design, \code{soptdesF}. 
}
\item{blk.N}{
integer, specifying  number of arrays \code{b} of resultant sequential optimal design, \code{soptdesF}.
}
\item{theta}{
numeric, representing  a function of the ratio of random array variance and random error variance. It takes any value between 0 and 1, inclusive. 
}
\item{soptdesF}{
matrix, a \code{2 x b} obtained sequential optimal or near-optimal block or row-column design.
}
\item{Optcrit}{
character specifying the optimality criteria to be used. Thus, \code{Optcrit} takes the letter \code{"A"}, \code{"MV"}, \code{"D"} and \code{"E"} for \code{A-}, \code{MV-}, \code{D-} and \code{E-}optimal or near-optimal block or row-column designs, respectively.
}
\item{strt}{
a non-negative integer, specifying number of added treatments/conditions to the initial design.  
}
\item{sary}{
a non-negative integer, specifying number of added arrays to the initial design. 
}
\item{dtype}{
character, specifying the design type. For block designs, \code{dtype = "blkd"} and for row-column deigns, \code{dtype = "rcd"}.
}
}
\details{
Detail discussions  concerning the constructions of graph can be found in \code{\link{igraph}} R package.
}
\value{
Returns the graphical layout of the resultant sequential optimal or near-optimal block or row-column design '\code{soptdesF}' on a separate pop-up window with the new edges (arrays) and vertices (treatments) added to the initial design  coloured in red and brown, respectively, for separation purpose. Furthermore, the function  \code{graphsoptd.mae} saves the graphical layout of the initial \code{des0} and resultant sequential optimal or near-optimal design in .pdf format in a working subdirectory. 
}
\author{
Dibaba Bayisa Gemechu, Legesse Kassa Debusho, and Linda Haines
}
\note{
When closing a pop-up window for graphical layout of the resultant designs (Graph plot),  if  the window is closed by 
clicking on the red button with "X" sign (top-right), the warning message  "Warning message: In \code{rm(list = cmd,} \code{envir = .tkplot.env):} object 'tkp \dots' not found"
will occur in R console irrespective of what command is executed next. To resolve this warning message, click 
on "\code{close}" menu that is located at the top-left of the graph plot pop-up window when closing this window. 
 }

\seealso{
\code{\link{igraph}}
}
\keyword{Graphical layout}
