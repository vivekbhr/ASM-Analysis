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
runDESeq(fcres,design, baseAllele = "CASTEiJ",topAllele = "129S1",fdrCutoff=0.01,
                     autodesigned = FALSE, SampleNames)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fcres}{
output of command "readFCres"
}
  \item{baseAllele}{
base level allele for comparison. eg. to count 129S1 over CASTEiJ, baseAllele = "CASTEiJ"
}

  \item{topAllele}{
top level allele for comparison. eg. to count 129S1 over CASTEiJ, topAllele = "129S1"
}

  \item{fdrCutoff}{
FDR cutoff for DESeq results
}

  \item{autodesigned}{
Logical (TRUE/FALSE) indicating if you used auto-designed sample matrix.
}

  \item{SampleNames}{
Data frame with sample names provided to readFCres command.
}


}
\details{

Returns a list object with DESeq results.

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
runDESeq(fcres,design, baseAllele = "CASTEiJ",topAllele = "129S1",fdrCutoff=0.01,
                     autodesigned = FALSE, SampleNames)
                          
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
