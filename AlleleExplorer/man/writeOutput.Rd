\name{writeOutput}
\alias{writeOutput}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  write back the ASexpression analysis with annotations.
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
  writeOutput(DESeqOutputList, fdrCutoff = 0.01, baseAllele = "CASTEiJ",topAllele = "129S1", 
                        autoAnnotate = TRUE, species = "Mus musculus", localGTF = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{DESeqOutputList}{
output of "runDESeq" command
}
  \item{fdrCutoff}{
FDR cutoff for DESeq results
}
    \item{baseAllele}{
base level allele for comparison. eg. to count 129S1 over CASTEiJ, baseAllele = "CASTEiJ"
}

  \item{topAllele}{
top level allele for comparison. eg. to count 129S1 over CASTEiJ, topAllele = "129S1"
}
 
    \item{autoAnnotate}{
TRUE if you want to fetch automatic annotations for output, FALSE if you want to provide your own GTF. (Default = TRUE)
}

  \item{species}{
The name of species to fetch annotation. Only required if using autoAnnotate. (Default = Mus musculus)
}

  \item{localGTF}{
Path of the local GTF file for annotation.
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
  
  writeOutput(myresultsList, fdrCutoff = 0.01, baseAllele = "CASTEiJ",topAllele = "129S1", 
                        autoAnnotate = TRUE, species = "Mus musculus", localGTF = NULL)
  
  # Writes back the files in present directory as [SampleName]_allelicBias_Output.txt
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
