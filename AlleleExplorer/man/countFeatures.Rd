\name{countFeatures}
\alias{countFeatures}

\title{
  perform featureCounts on the given set of samples.
}
\description{

}
\usage{
countFeatures(SampleInfo,AutoParam="RNASeq",AnnotationFile="~/Desktop/myAnnotation.GTF",nthreads=10,outFileName=NULL,...)
}

\arguments{
  
  \item{SampleInfo}{
    information of Samples in a data.frame with four columns, containing SampleID(filename with Path), a Name (eg, TF1), 
    replicate number (1,2,3..) and Allele (mat/pat).
  }
  \item{AutoParam}{
    Could be "RNASeq","ChIPSeq" or NULL, to specify your own parameters to featureCounts.
  }
  \item{AnnotationFile}{
    Path of GTF or SAF file for annotation. In case of AutoParam="RNASeq", GTF input is required, 
    while in "ChIPSeq" a SAF input is needed. In case of AutoParam=NULL, you need to specify your own format, see
    featureCounts function for details.
  }
  
  \item{nthreads}{
    number of threads to compute.
  }
  
  \item{outFileName}{
    Filename to write your output. Null for capturing output as a list.
  }
  
  \item{...}{
    Any other parameter to featureCounts. See Rsubread featureCounts help file for details.
  }
  
  
  
}
\details{

}
\value{

}
\references{

}
\author{
  Vivek Bhardwaj
}
\note{
}

\examples{
  
  1.. Make a SampleInfo file.
  
  Samples <- data.frame(samp=c("/path/to/sample/sample1.txt","/path/to/sample/sample1.txt",
                               "/path/to/sample/sample1.txt","/path/to/sample/sample1.txt"),
                        name=c("TF1","TF1","TF1","TF1"), rep=c(1,1,2,2), 
                        allele= c("mat","pat","mat","pat"))
  Samples
  samp name rep allele
  1 /path/to/sample/sample1.txt  TF1   1    mat
  2 /path/to/sample/sample1.txt  TF1   1    pat
  3 /path/to/sample/sample1.txt  TF1   2    mat
  4 /path/to/sample/sample1.txt  TF1   2    pat
  

  2.. Run the function
  countFeatures(SampleInfo=Samples,AutoParam="RNASeq",AnnotationFile="~/Desktop/myAnnotation.GTF",nthreads=10,outFileName=NULL,...)
  
  
}

