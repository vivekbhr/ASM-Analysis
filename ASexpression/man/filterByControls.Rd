\name{filterByControls}
\alias{filterByControls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
filter out the genes which are also biased in control samples.
}
\description{
Assuming the genes biased in controls are not of interest, this script filters them out.
}
\usage{
filterByControls(DESeqResultList, matchFile, padjCutoff = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DESeqResultList}{
output of "pullResults" command.
}
  \item{matchFile}{
a control to test matching file (see example).
}
  \item{padjCutoff}{
diff-exp cutoff for exclusion.
}
}
\details{
output is a list with data frame for each TF. can be written on a file.
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
matchfile = data.frame(test = c("MOF","Msl1","NSL3"), controls = c("scr","Scr","scr"))

diff.unique <- filterByControls(ase.results,matchfile)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
