\name{pathwayEnrichment}
\alias{pathwayEnrichment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot pathway enrichment of the DESeq output
}
\description{
takes as input the output of "runDESeq" command and gives enriched pathway for each TF as a list object.
}
\usage{
pathwayEnrichment(DESeqOutputList,organism="mmu",fdrCutoff=0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DESeqOutputList}{
output of "runDESeq" command
}
  \item{organism}{
three word kegg code for organism, eg. mmu, hsa. default = "mmu"
}
  \item{fdrCutoff}{
fdrCutoff for DESeq outputs to plot the results
}
}
\details{
Output is just a list of enriched pathway for each TF in your sample. Later can be written in a file.
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
pathways <- pathwayEnrichment(DESeqOutputList,organism="mmu",fdrCutoff=0.01)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
