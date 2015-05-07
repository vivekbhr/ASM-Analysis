\name{readFCres}
\alias{readFCres}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
read featurecount output with Allele-specific reads
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
readFCres(CountFile, ContSampleNames, TestSampleNames)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CountFile}{
Raw featureCounts output file, tab seperated with headers.
}
  \item{ContSampleNames}{
A text file with the name of controls(scrambled/Wildtype) samples, seperated by newline.
}
  \item{TestSampleNames}{
A text file with the name of Test(knockdown) samples, seperated by newline.
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
fcres.ase <- readFCres(CountFile = "~/Desktop/featCount_genes_genesCollapsed_RNAseqSamples.txt",ContSampleNames = "~/Desktop/contnames.txt",TestSampleNames = "~/Desktop/samplenames.txt")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
