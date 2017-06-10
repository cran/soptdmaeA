\name{soptdmaeA-internal}
\alias{soptdmaeA-internal}
\alias{seqAoptbrcd.maeA}
\alias{seqMVoptbrcd.maeA}
\alias{seqDoptbrcd.maeA}
\alias{seqEoptbrcd.maeA}
%\alias{desoI}
\title{
Internal functions
}
\description{
Functions for internal usage only.
}
\usage{
 
## Computes A-optimal or near-optimal block or row-column designs
## using array exchange algorithm
seqAoptbrcd.maeA(trt.N, blk.N, theta, nrep, strt, sary, des0, dtype)
 
## Computes MV-optimal or near-optimal block or row-column designs
## using array exchange algorithm
seqMVoptbrcd.maeA(trt.N, blk.N, theta, nrep, strt, sary, des0, dtype)
 
## Computes A-optimal or near-optimal block or row-column designs
## using array exchange algorithm
seqDoptbrcd.maeA(trt.N, blk.N, theta, nrep, strt, sary, des0, dtype)
 
## Computes A-optimal or near-optimal block or row-column designs
## using array exchange algorithm
seqEoptbrcd.maeA(trt.N, blk.N, theta, nrep, strt, sary, des0, dtype)
% 
%## Computes A-optimal or near-optimal block or row-column designs
%## using array exchange algorithm
%##desoI(des0)


 
}
\arguments{
  \item{trt.N}{
integer, specifying number of treatments \code{v} of initial design, \code{des0}. 
}
  \item{blk.N}{
integer, specifying number of arrays \code{b} of initial design, \code{des0}.
}
  \item{theta}{
numeric, representing  a function of the ratio of random array variance and random error variance. It takes any value between 0 and 1, inclusive. 
}
  \item{nrep}{
integer, specifying number of replications of the optimization procedure. 
}
  \item{strt}{
a non-negative integer, specifying number of added treatments/conditions to the initial design.  
}
  \item{sary}{
a non-negative integer, specifying number of added arrays to the initial design. 
}
  \item{des0}{
matrix, a \code{2 x blk.N} or \code{blk.N x 2} initial block or row-column design. The initial design must be treatment connected and the number of treatments and arrays should also coincides with \code{trt.N} and \code{blk.N} inserted by the user, if this conditions are not satisfied, the package will stop running with an error message. 
}
  \item{dtype}{
character, specifying the design type. For block designs, \code{dtype = "blkd"} and for row-column deigns, \code{dtype = "rcd"}.
}
}
\details{
These functions are handled via a generic function \code{\link{soptdmaeA}}. Please refer to the \code{\link{soptdmaeA}} documentation for details. 
}

\references{
Debusho, L. K., Gemechu, D. B., and Haines, L. M. (2016).  Algorithmic construction of optimal block designs for two-colour cDNA microarray experiments using the linear mixed model. Under review.

Gemechu D. B., Debusho L. K. and Haines L. M. (2014). A-optimal designs for two-colour cDNA microarray experiments using the linear mixed effects model. \emph{Peer-reviewed Proceedings of the Annual Conference of the South African Statistical Association for 2014 (SASA 2014), Rhodes University, Grahamstown, South Africa}. pp 33-40, ISBN: 978-1-86822-659-7.
}
\author{
Dibaba Bayisa Gemechu, Legesse Kassa Debusho, and Linda Haines
}
\seealso{
\code{\link{soptdmaeA}}
}
\keyword{Sequential A-optimal block designs}
\keyword{Sequential D-optimal block designs}
\keyword{Sequential E-optimal block designs}
\keyword{Sequential MV-optimal block designs}
\keyword{Sequential A-optimal row-column designs}
\keyword{Sequential D-optimal row-column designs}
\keyword{Sequential E-optimal row-column designs}
\keyword{Sequential MV-optimal row-column designs}
\keyword{Microarray experiment} 
\keyword{Array exchange algorithm} 
