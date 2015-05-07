\name{runDESeq}
\alias{runDESeq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
run DESeq on the output of readFCres.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
runDESeq(fcres, design.cont, design.test)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fcres}{
output of command "readFCres"
}
  \item{design.cont}{
control design file with two columns (see example).
}
  \item{design.test}{
test design file with two columns (see example).
}
}
\details{

The design matrix should have row.names = colnames(countOutput), and two columns:
sample = name of sample
allele = name of allele (maternal/paternal)

The column names in readFCres output should correspond to the description in the design matrix in the same order.
For example, if colnames(fcres) = "control1_maternal", "control2_maternal","control1_paternal", "control2_paternal" ; then design matrix should have sample = c("cont1", "cont1", "cont2","cont2") ; and allele = c("maternal","paternal","maternal","paternal").


returns a list object with summarizedExp class objects.

}
\value{
%%
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Vivek Bhardwaj
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
testdesign <- data.frame(row.names = colnames(fcres.ase$test),
                          sample = rep(rep(c("MOF","Msl1","NSL3"),each = 3),2),
                          allele = rep(c("129S1","CASTEiJ"), each = 9))
contdesign <- data.frame(row.names = colnames(fcres.ase$cont),
                          sample = rep(c("scr","Scr"),6), 
                          allele = rep(c("129S1","CASTEiJ"),each=6))
runDESeq(fcres="readFCresoutput", design.cont="contdesign", design.test="testdesign")
                          
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
