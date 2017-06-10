\name{fixparsoptd.mae}
\alias{fixparsoptd.mae}
\title{
Sets parametric values and initial design
}
\description{
Creates a \code{\link{tcltk}} window that allow the users to set or fix values for the parametric combinations and initial optimal or near-optimal design to compute sequential optimal or near-optimal block or row-column designs.}
\usage{
fixparsoptd.mae(Optcrit)
}
\arguments{
  \item{Optcrit}{
character, specifying the optimality criteria to be used. \code{Optcrit} takes the letter \code{"A"}, \code{"MV"}, \code{"D"} and \code{"E"} for \code{A-}, \code{MV-}, \code{D-} and \code{E-}optimal or near-optimal designs, respectively.
}
}
\details{
\code{fixparsoptd.mae} creates a pop-up tcltk window that allow the users to set 
the parametric combinations  and initial optimal or near-optimal designs to compute optimal or near-optimal block or row-column designs. The parameters include the number of treatments \code{trt.N}, number arrays  \code{blk.N}, theta value  \code{theta}, number of replications of the optimization procedure \code{nrep}, number of added treatments \code{strt}, number of added arrays \code{sary} and a button that allows users to insert initial optimal or near-optimal design \code{des0}  by typing into a pop-up sheet that will occur when clicking on the 'Insert' button. The users are not restricted to a specific dimension of the initial design, thus,  the dimension of the initial design \code{des0} can be \code{2 x blk.N} or \code{blk.N x 2}. The default initial design \code{des0} is a \code{3 x 2} loop/cyclic design. 

Furthermore, on this window, the checkbox options that allow the users to choose whether or not they need to have the graphical 
layout of the resultant sequential optimal or near-optimal block or row-column design, to make a choice between the two design types (block design and row-column designs) and to print the summary of the resultant sequential optimal or near-optimal design on R-console directly are available. 

After setting all the required parametric combinations and selecting the algorithm of interest, 
clicking on the search button on the set parametric combinations tcltk window, similar to the results that 
can be obtained when using the function \code{\link{soptdmaeA}}, the summary of the resultant sequential optimal or near-optimal design will be saved in the current working directory in .csv format. The graphical layout of resultant sequential optimal or near-optimal design together with that of the initial optimal or near-optimal design  will also be displayed on R console.
}
\value{
The \code{fixparsoptd.mae} function creates a pop-up tcltk window that allow the users to set 
the parametric combinations to compute sequential optimal or near-optimal block or row-column designs.
}
\author{
Dibaba Bayisa Gemechu, Legesse Kassa Debusho, and Linda Haines
}
 
\seealso{
\code{\link{soptdmaeA}}, \code{\link{mmenusoptd.mae}}, \code{\link{tcltk}}, \code{\link{TkWidgets}}
}
 