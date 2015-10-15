\name{createDesignMatrix }
\alias{createDesignMatrix }
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
automatically create design matrix from the countgenes output.
}
\description{
%
}
\usage{
createDesignMatrix(fcres, SampleNames, baseAllele = "CASTEiJ")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fcres}{
output of readFCres command
}
  \item{SampleNames}{
Data frame with sample names provided to readFCres command.
}
  \item{baseAllele}{
base level allele for comparison. eg. to count 129S1 over CASTEiJ, baseAllele = "CASTEiJ"
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
# make SampleNames data frame
# first run readFCres command
fcres <- readFCres(CountFile,SampleNames)
# then create design matrix
createDesignMatrix <- function(fcres, SampleNames, baseAllele = "CASTEiJ")


## you can also create design matrix on your own, it has to be either a data.frame or a list with names = name of sample (same as col1 in SampleNames). Each element of the list looks like this

row.names       condition  allele 
Msl1_1_129S1          KD   129S1  
Msl1_2_129S1          KD   129S1  
Msl1_3_129S1          KD   129S1  
Msl1_1_CASTEiJ        KD CASTEiJ  
Msl1_2_CASTEiJ        KD CASTEiJ  
Msl1_3_CASTEiJ        KD CASTEiJ  
ScrM_1_129S1     Control   129S1  
ScrM_2_129S1     Control   129S1  
ScrM_3_129S1     Control   129S1  
ScrM_1_CASTEiJ   Control CASTEiJ  
ScrM_2_CASTEiJ   Control CASTEiJ  
ScrM_3_CASTEiJ   Control CASTEiJ  


with row.names = colnames(fcres$sample) and two columns with compared samples and alleles. 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
