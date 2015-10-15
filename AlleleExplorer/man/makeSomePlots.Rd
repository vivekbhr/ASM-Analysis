\name{makeSomePlots}
\alias{makeSomePlots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot PCA of samples and the number of ASE genes after filtering
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makeSomePlots(DESeqOutputList,baseAllele = "CASTEiJ",topAllele = "129S1",fdrCutoff=0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DESeqOutputList}{
output of "runDESeq " command
}
  \item{baseAllele}{
base level allele for comparison. eg. to count 129S1 over CASTEiJ, baseAllele = "CASTEiJ"
}

  \item{topAllele}{
top level allele for comparison. eg. to count 129S1 over CASTEiJ, baseAllele = "129S1"
}

  \item{fdrCutoff}{
FDR cutoff for DESeq results
}

}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
diffExpPlots <- makeSomePlots(DESeqOutputList,baseAllele = "CASTEiJ",topAllele = "129S1",fdrCutoff=0.01)

# diffexpPlots shows the number of genes detected with bias for each allele. On top of that you get PCA plot on screen
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
