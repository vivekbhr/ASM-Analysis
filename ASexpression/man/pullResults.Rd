\name{pullResults}
\alias{pullResults}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
pull the results from DESeqOutput list.
}
\description{
takes as input the output of runDESeq command, and fetches the results using the given compare group.
}
\usage{
pullResults(DESeqOutputList, maternal = "129S1", paternal = "CASTEiJ", compareGroup = "MatOverPat", design.cont, design.test)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DESeqOutputList}{
output of runDESeq command
}
  \item{maternal}{
name of maternal allele.
}
  \item{paternal}{
name of paternal allele.
}
  \item{compareGroup}{
how to compare. One of "MatOverPat" or "PatOverPat"
}
  \item{design.cont}{
control design matrix. same as the one used for runDESeq command.
}
  \item{design.test}{
test design matrix. same as the one used for runDESeq command.
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
ase.results <- pullResults(DESeqOutputList = ase.deseq,maternal = "129S1",paternal = "CASTEiJ",
                          compareGroup = "MatOverPat",
                           design.cont = contdesign, design.test = testdesign)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
