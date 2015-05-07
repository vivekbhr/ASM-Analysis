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
makeSomePlots(DESeqOutputList, ControlFilteredOutput, compareGroup = "MatOverPat")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DESeqOutputList}{
output of "runDESeq " command
}
  \item{ControlFilteredOutput}{
output of filterbyControls command.
}
  \item{compareGroup}{
same as the one used for "pullResults" command
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
makeSomePlots(DESeqOutputList = ase.deseq,ControlFilteredOutput = diff.unique,compareGroup = "MatOverPat")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
