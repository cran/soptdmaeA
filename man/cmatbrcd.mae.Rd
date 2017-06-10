\name{cmatbrcd.mae}
\alias{cmatbrcd.mae}
\title{
Computes the treatment information matrix 
}
\description{
Computes the information matrix (C-matrix) for treatment effects under either the linear fixed effects model  or the linear mixed effects model setting for a given block or row-column design.
}
\usage{
cmatbrcd.mae(trt.N, blk.N, theta, des, dtype)
}
\arguments{
  \item{trt.N}{
integer, specifying number of treatments \code{v} of design \code{des}. 
}
  \item{blk.N}{
integer, specifying number of arrays (blocks or columns) \code{b} of design \code{des}.
}
  \item{theta}{
numeric, representing  a function of the ratio of random array variance and random error variance. It takes any value between 0 and 1, inclusive. 
}
  \item{des}{
matrix, a \code{2 x b} block or row-column design with \code{b} arrays and \code{v} treatments.}
  \item{dtype}{
character, specifying the design type. For block designs, \code{dtype = "blkd"} and for row-column deigns, \code{dtype = "rcd"}.
}
}
\value{
Returns a \code{v x v} treatment information matrix (C-matrix).}
\references{
Debusho, L. K., Gemechu, D. B., and Haines, L. M. (2016).  Algorithmic construction of optimal block designs for two-colour cDNA microarray experiments using the linear mixed model. Under review.

Gemechu D. B., Debusho L. K. and Haines L. M. (2014). A-optimal designs for two-colour cDNA microarray experiments using the linear mixed effects model. \emph{Peer-reviewed Proceedings of the Annual Conference of the South African Statistical Association for 2014 (SASA 2014), Rhodes University, Grahamstown, South Africa}. pp 33-40, ISBN: 978-1-86822-659-7.
}
\author{
Dibaba Bayisa Gemechu, Legesse Kassa Debusho, and Linda Haines
}

\seealso{
\code{\link{soptdmaeA}}, \code{\link{fixparsoptd.mae}}
}
\examples{

##Information matrix

     trt.N <- 4 
     blk.N <- 4 
     theta <- 0.3 
     dsgn <- rbind(1:4,c(2:4,1))
     dtype <- "rcd"

     cmatbrcd.mae(trt.N = 4, blk.N = 4, theta = 0.2, des = dsgn, dtype = "rcd")
}
\keyword{Information matrix}
\keyword{C-matrix} 
