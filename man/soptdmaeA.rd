\name{soptdmaeA}
\alias{soptdmaeA}
\alias{soptdmaeA.default}
\alias{print.soptdmaeA}
\alias{summary.soptdmaeA}
\alias{print.summary.soptdmaeA}
\title{
Sequential optimal designs for two-colour cDNA microarray experiments
}
\description{
Used to compute sequential A-, MV-, D- or E-optimal or near-optimal block and row-column designs for two-colour cDNA microarray experiments under either the linear fixed effects model or the linear mixed effects model settings using the array exchange algorithms of Debusho, Gemechu and Haines (2016).}
\usage{
soptdmaeA(trt.N, blk.N, theta, nrep, strt, sary, des0, dtype, Optcrit = "", ...)


\method{soptdmaeA}{default}(trt.N, blk.N, theta, nrep, strt, sary, des0, dtype, Optcrit = "",...)
\method{print}{soptdmaeA}(x, ...)
\method{summary}{soptdmaeA}(object, ...)
}
\arguments{
  \item{trt.N}{
integer, specifying number of treatments \code{v} of initial design, \code{des0}. 
}
  \item{blk.N}{
integer, specifying number of arrays (blocks or columns) \code{b} of initial design, \code{des0}.
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
character, specifying the design type. For block designs, \code{dtype = "blkd"} and for row-column deigns, \code{dtype = "rcd"}.}
  \item{Optcrit}{
character, specifying the optimality criteria to be used. \code{Optcrit} takes the letter \code{"A"}, \code{"MV"}, \code{"D"} and \code{"E"} for \code{A-}, \code{MV-}, \code{D-} and \code{E-}optimal or near-optimal block or row-column designs, respectively.
}
  \item{x}{
the object to be printed.
}
  \item{object}{
an object of class \code{"soptdmaeA"}.
}
  \item{\dots}{
not used.
}
}
\details{
\code{soptdmaeA} computes sequential optimal or near-optimal block or row-column designs for the two-colour cDNA microarray experiments  where the interest is in a comparison of all possible elementary treatment contrasts for a given initial optimal or near-optimal designs. The function computes sequential A-, MV-, D- and E-optimal or near optimal block or row-column designs via calling of four sub-functions \code{\link{seqAoptbrcd.maeA}}, \code{\link{seqMVoptbrcd.maeA}}, \code{\link{seqDoptbrcd.maeA}}, and \code{\link{seqEoptbrcd.maeA}}, respectively. These functions uses the array exchange algorithm of Debusho, Gemechu and Haines (2016). Thus, once the parametric combinations of interest are sated, these functions will first compute, randomly, a new connected initial design with a new number of arrays and, optionally, a new number of treatments. Then they perform the array exchange procedure through deletion and addition of candidate arrays at a time and selects a design with best array exchange with respect to the optimality criterion value. The candidate arrays are lists of possible arrays with different treatment combinations and their lists are dependent of the number of arrays and treatments added to the initial optimal or near-optimal design. For example, if only one treatment and one array are to be added to the initial optimal or near-optimal design, then the candidate arrays will be only those arrays that consists of a new treatment together with the old treatments in the initial optimal or near-optimal design with or without considering their position within the array for row-column or block designs, respectively.    

The minimum value of \code{trt.N} and \code{blk.N} is 3 and \code{trt.N} should be less than or equal to \code{blk.N - 1}. Thus, the least initial design should be of a design with 3 number of treatments and number of arrays. The minimum number of \code{sary} and \code{strt} are 1 and 0, respectively, and \code{sary} should be greater than or equal to \code{strt}. 
The linear fixed effects model results for given parametric combinations and initial design are obtained by setting \code{theta = 0.0}.

\code{nrep} takes a value of greater than or equal to 1. However, to ensure optimality of the resultant design, for \code{sary - strt > 0},  
the \code{nrep} should be greater than or equal to 10. In addition, as \code{trt.N} or \code{blk.N} or \code{sary} and/or \code{strt} or all of them increase, 
to ensure optimality of resultant design, it is advised to further increase the value of \code{nrep}
up to greater than or equal to 50. However, it has to be noted that as \code{trt.N} or \code{blk.N} or
 \code{nrep} or all of them increase, computer time required to generate sequential optimal or near-optimal design increases.

}
\value{
Returns the initial and resultant sequential A-, MV-, D- or E-optimal or near-optimal block or row-column design with their corresponding score value and parametric combination 
saved in excel file in a working directory. In addition, the function \code{soptdmaeA} displays the graphical layout of the initial and resultant 
optimal or near-optimal block or row-column designs. Specifically: 

\item{call}{the method call.}         
\item{v}{number of treatments of obtained sequential design.}
\item{b}{number of arrays  of obtained sequential design.}
\item{theta}{theta value.}
\item{nrep}{number of replications of the optimization procedure.}  
\item{strt}{number of added treatments.}                          
\item{sary}{number of added arrays.}                          
\item{Optcrit}{optimality criteria.} 
\item{optdes0}{a \code{2 x blk.N} initial optimal or near-optimal block or row-column design.}
\item{optcrtsv0}{score value of the optimality criteria \code{'Optcrit'} of the initial optimal or near-optimal block or row-column design \code{'optdes0'}.}
\item{soptdesF}{a \code{2 x blk.N} obtained sequential optimal or near-optimal block or row-column design.}
\item{soptcrtsv}{score value of the optimality criteria \code{'Optcrit'} of the resultant sequential optimal or near-optimal block design \code{'soptdesF'}.}
\item{file_loc, file_loc2}{location where the summary of the resultant optimal or near-optimal block design is saved in .csv format.}
\item{equireplicate0}{logical value indicating whether the initial optimal or near-optimal block or row-column design is equireplicate or not.}
\item{vtrtrep0}{vector of treatment replication of the initial optimal or near-optimal block or row-column design.}     
\item{equireplicate}{logical value indicating whether the resultant sequential optimal or near-optimal block or row-column design is equireplicate or not.}
\item{vtrtrep}{vector of treatment replication of the resultant sequential optimal or near-optimal block or row-column design.}     
\item{Cmat}{the C-matrix or  treatment information matrix of the  obtained sequential optimal or near-optimal block or row-column design.} 

The output also includes graphical layouts of the initial and resultant sequential optimal or near-optimal block or row-column design.   The new edges (arrays) and vertices (treatments) added to the initial design are coloured in red and brown, respectively, for identification purpose. 

NB: The function \code{soptdmaeA} also saves the summary of the initial and resultant sequential optimal or near-optimal block or row-column design in .csv format in the working directory. 
Furthermore, the function reports only one final sequential optimal or near-optimal block or row-column design, however, there is a possibility 
of more than one sequential optimal or near-optimal block or row-column designs for a given parametric combination. 
The function \code{\link{graphsoptd.mae}} can be used to view and rearrange the graphical layout of the resultant 
sequential optimal or near-optimal block or row-column design on \code{tcltk} window. Alternative to the function \code{soptdmaeA}, a
GUI tcltk window can be used to generate sequential optimal or near-optimal block or row-column designs, see \code{\link{mmenusoptd.mae}} and \code{\link{fixparsoptd.mae}}.   

}
\references{
Debusho, L. K., Gemechu, D. B., and Haines, L. M. (2016).  Algorithmic construction of optimal block designs for two-colour cDNA microarray experiments using the linear mixed model. Under review.

Gemechu D. B., Debusho L. K. and Haines L. M. (2014). A-optimal designs for two-colour cDNA microarray experiments using the linear mixed effects model. \emph{Peer-reviewed Proceedings of the Annual Conference of the South African Statistical Association for 2014 (SASA 2014), Rhodes University, Grahamstown, South Africa}. pp 33-40, ISBN: 978-1-86822-659-7.
}
\author{
Dibaba Bayisa Gemechu, Legesse Kassa Debusho, and Linda Haines
}
\seealso{
\code{\link{mmenusoptd.mae}}, \code{\link{fixparsoptd.mae}}
}

\examples{
  \donttest{
  ##To obtain sequential A-optimal or near-optimal block design for a given
  ##initial A-optimal or near-optimal block design, set
  
  trt.N <- 3  #Number of treatments
  blk.N <- 3  #Number of blocks
  theta <- 0  #theta value
  nrep <- 10  #Number of replications
  strt <- 2   #Number of added treatments
  sary <- 3   #Number of added arrays
  des0 <- rbind(1:3, c(2, 3, 1)) #Initial design
  dtype = "blkd" #Design type
  Optcrit <- "A"   #Optimality criteria

  seqAoptbd <- soptdmaeA(trt.N = 3, blk.N = 3, theta = 0, nrep = 10, 
                            strt = 2, sary = 3, des0, dtype = "blkd", Optcrit = "A")
  
  summary(seqAoptbd)

  ##To obtain sequential A-optimal or near-optimal row-column design for a given
  ##initial A-optimal or near-optimal row-column design des0 (stated above), set
  
  dtype = "rcd" #Design type

  seqAoptrcd <- soptdmaeA(trt.N = 3, blk.N = 3, theta = 0, nrep = 10, 
                            strt = 2, sary = 3, des0, dtype = "rcd", Optcrit = "A")
  
  summary(seqAoptrcd)
}
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
