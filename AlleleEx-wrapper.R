#!/usr/bin/env Rscript
## A wrapper for all R functions of AlleleExplorer
## Copyright 2015 Vivek Bhardwaj (bhardwaj@ie-freiburg.mpg.de). Licence: GPLv3.

options(warn = -1)
suppressPackageStartupMessages({
  require(AlleleExplorer)
  require(argparser)
})


## Add arguments
p <- arg_parser("A wrapper for all R functions of AlleleExplorer")
p <- add_argument(p,"--sampleSheet",help = "The sample sheet to parse")
p <- add_argument(p,"--refAllele",help = "Reference allele (to count fold changes/enrichment on) ")
p <- add_argument(p,"--fdr",help = "FDR cutoff")
p <- add_argument(p,"--outdir",help = "Output directory")
p <- add_argument(p,"--gtf",help = "GTF annotation file for RNA-Seq")
p <- add_argument(p,"--exclude",help = "Any chromosome to exclude from output")
p <- add_argument(p,"--threads",help = "Number of cores to use (leave empty for one)")
p <- add_argument(p,"--qcplots",flag = TRUE, help = "Whether you want to make QC plots")
## Parse the arguments
argv <- parse_args(p)

sheet <- argv$sampleSheet
ref <- argv$refAllele
mode <- argv$mode
fdr <- argv$fdr
out <- argv$outdir
gtf <- argv$gtf
exclude <- argv$exclude
threads <-argv$threads
qcplots <-argv$qcplots

### ----------------------------------  Write and run the functions ---------------------------

samp <- read.csv(sheet,header=TRUE,
                 colClasses = c("character","character","factor","factor","factor","character"))

## ChIP samples
if(any(grepl("chip",samp[,1]))){
  message("ChIP-Seq samples found. Running allele-specific binding analysis")
  chipsamp <- dplyr::filter(samp,expType == "chip",sampleType == "test")
  tflist <- as.character(unique(chipsamp$tf))

  chipCountObject <- readfiles_chip(csvFile = sheet, refAllele = ref)
  if(qcplots){
    message("making QC plots")
    makeQCplots_chip(chipCountObject,out)
  } else message("skipping QC plots")
  message("filtering chip by input")
  chipCountObject <- filterByInput_chip(chipCountObject)
  message("normalizing windows by tmm")
  chipCountObject <- tmmNormalize_chip(chipCountObject, plotfile = paste0(out,"/TMM_normalizedCounts.pdf"))
  for(tf in tflist){
    tryCatch({
      message(paste0("extracting Db regions for : ",tf))
      chipResultObject <- getDBregions_chip(chipCountObject,
                                            plotfile = paste0(out,"/Dispersions.pdf"), tfname = tf)
      message(paste0("Writing output for : ",tf))
      writeOutput_chip(chipResultObject, outfileName = paste0(out,"/ASE-result_",tf) , fdr = fdr, annotation = FALSE,
                       Txdb = TxDb.Mmusculus.UCSC.mm9.knownGene, Orgdb = org.Mm.eg.db)

      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
        )
      }
} else {
  message("ChIP-Seq samples not found.")
}

## RNA samples
if(any(grepl("rna",samp[,1]))){
  message("RNA-Seq samples found. Running allele-specific expression analysis")
  rnasamp <- dplyr::filter(samp,expType == "rna",sampleType == "test")
  tflist <- as.character(unique(rnasamp$tf))
  message("Counting features")
  rnaCountObject <- countFeatures_rna(csvFile = sheet ,AnnotationFile = gtf, nthreads = threads,
                                      outfileName = paste0(out,"/ASE_counts.tab"),
                                      refAllele = ref,userParam = FALSE)
  for(tf in tflist){
        tryCatch({
          message(paste0("running differential expression for : ",tf))
          rnaResultObject <- alleleDiff_rna(rnaCountObject,fdrCutoff = fdr,tfname = tf)
          plotResults_rna(rnaResultObject,outfile = paste0(out,"resultPlots_RNA_",tf,".pdf"),
                          barplot = TRUE, excludeChr = exclude)
          message(paste0("writing output for : ",tf))
          writeOutput_rna(rnaResultObject,annotateFrom = "dataset",excludeChr = exclude,
                          fdrCutoff = fdr,outfileName = paste0(out,"/ASE-result_",tf) )
          }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
          )
        }
  } else {
    message("RNA-Seq samples not found.")
}

message("DONE..!!")
