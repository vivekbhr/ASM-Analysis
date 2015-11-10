
# Functions to perform allele-sp differential expression
# (c) Vivek Bhardwaj (bhardwaj@ie-freiburg.mpg.de)
# First created: 13th May 2015

#' Read the Files and Count features for the RNA-Seq Samples
#'
#' @param csvFile csvfile with sample information
#' @param AnnotationFile GTF file
#' @param nthreads Number of threads
#' @param outfileName Specify output file
#' @param refAllele Reference Allele (name must match the name in samplesheet)
#' @param userParam Make it TRUE if you want to enter your own featurecounts parameters
#' @return rnaCountObject : a list with counts and sampleinfo
#' @examples
#' countFeatures_rna(csvFile = "testBAMs/testSampleSheet.csv",AnnotationFile,nthreads,
#' outfileName=NULL,refAllele = "pat",userParam = FALSE,...)
#' 

countFeatures_rna <- function(csvFile = "testBAMs/testSampleSheet.csv",AnnotationFile,nthreads,
                              outfileName=NULL,refAllele = "pat",userParam = FALSE,...){
  
  # Parse Sample sheet
  samp <- read.csv(csvFile,header=TRUE,
                   colClasses = c("character","character","factor","factor","factor","character"))
  if(any(grepl("rna",samp[,1]))){
    message("Extracting sample Info")
    samp <- samp[which(grepl("rna",samp[,1])),]
  } else stop("No RNA Samples found for diffBinding test. Check sample sheet")
  
  ## make design matrix for DESeq2
  design <- data.frame(row.names = samp[,2] , samp[,3:5]) ## Add: A warning if rownames not uniq
  colnames(design) <- c("condition","allele","tf")
  design$allele <- relevel(design$allele,refAllele)
  design$condition <- relevel(design$condition,"control")
  
  # Set the TF corresponding to control as the reference
  reftf <- as.character(unique(design[which(design$condition == "control"),"tf"]))
  if(length(reftf) != 1) stop("Condition 'control' has multiple tfs associated. 
                              Make sure you have only one tf as control")
  design$tf <- relevel(design$tf,reftf)
  
  # finally save the name of ref and alt allele
  altAllele <- unique(design[which(design$allele != refAllele),"allele"])
  if(length(altAllele) != 1) stop("More than one alternative alleles in samplesheet?")
  alleleinfo <- data.frame(refAllele = refAllele, altAllele = altAllele)
  
  ## Count the reads
  bam.files <- samp[,6]
  message("Using Default options for Allele-Specific Expression")
  if(userParam){
    message("Counting reads using user-provided parameters")
    rnaCountObject <- Rsubread::featureCounts(files = bam.files,annot.ext = AnnotationFile,
                                              isGTFAnnotationFile=TRUE,...)
  } else {
    message("Counting reads using standard parameters")
    rnaCountObject <- Rsubread::featureCounts(files = bam.files,annot.ext = AnnotationFile,
                                              isGTFAnnotationFile=TRUE,useMetaFeatures=TRUE,
                                              isPairedEnd=TRUE,requireBothEndsMapped=FALSE,checkFragLength=FALSE,
                                              nthreads=nthreads,strandSpecific=2,minMQS=0,countPrimaryAlignmentsOnly=TRUE)
  }
  # Save Fcount output 
  colnames(rnaCountObject$counts) <- samp[,2] # change colnames to samplenames
  rnaCountObject$targets <- samp[,2] # changed targets(just another df) to samplenames
  if(!(is.null(outfileName))){
    write.table(rnaCountObject$counts,outfileName,sep="\t",row.names=F,quote=F)
  }
  
  ## Make the output
  rnaCountObject$design <- design
  rnaCountObject$alleleinfo <- alleleinfo
  return(rnaCountObject)
}


#' Run DESeq2 with interaction design
#'
#' @param rnaCountObject output from countFeatures_rna script
#' @param fdrCutoff FDR cutoff for differential expression
#' @param tfname Name of TF to get the differential expression (must match the name in samplesheet)
#' @return rnaResultObject
#' @examples
#' alleleDiff_rna(rnaCountObject,fdrCutoff = 0.01,tfname = "mof")
#' 

alleleDiff_rna <- function(rnaCountObject,fdrCutoff = 0.01,tfname = "mof"){ 
  
  countdata <- rnaCountObject$counts
  designmat <- rnaCountObject$design
  alt <- as.character(rnaCountObject$alleleinfo$altAllele)
  # Run DESeq (all info already in the rnaCountObject)
  dds <- DESeq2::DESeqDataSetFromMatrix(countdata,designmat,
                                        design = ~ tf * allele)
  dds <- DESeq2::DESeq(dds,betaPrior = FALSE)  
  # Get results
  message("DESeq finished. Now extracting results.")
  
  query <- paste0("tf",tfname,".allele",alt)
  ddr = DESeq2::results(dds, name = query, cooksCutoff = FALSE,alpha = fdrCutoff)
  
  # Return outputs
  rnaCountObject$DEdataSet <- dds
  rnaCountObject$DEresult <- ddr
  return(rnaCountObject) ## Now it's called "rnaResultObject"
}

#' Make QC and result plots for DESeq2 output
#'
#' @param rnaResultObject Output of alleleDiff_rna
#' @param outfile Output file to write back the result
#' @param barplot TRUE/FALSE whether you want the barplot for allele-biased gene numbers
#' @param excludeChr The chromosome to exclude from the barplot (if barplot = TRUE)
#' @return A pdf file with QC and results plots
#' @examples
#' plotResults_rna(rnaResultObject,outfile = "resultPlots_RNA.pdf")
#' 
## Function requires : vsn, DESeq2 > 1.10, pheatmap, RColorBrewer, ggplot2

plotResults_rna <- function(rnaResultObject,outfile = "resultPlots_RNA.pdf", barplot = TRUE,
                            fdrCutoff = 0.01, excludeChr = "chr12" ){
  dds <- rnaResultObject$DEdataSet
  ddr <- rnaResultObject$DEresult
  rld <- DESeq2::rlog(dds) #rlog transform
  vsd <- DESeq2::varianceStabilizingTransformation(dds) # vst transform
  
  pdf(outfile)
  # Plot transformed counts
  notAllZero <- (rowSums(counts(dds))>0)
  par(mfrow = c(1,3))
  vsn::meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1))
  vsn::meanSdPlot(assay(rld[notAllZero,]))
  vsn::meanSdPlot(assay(vsd[notAllZero,]))
  par(mfrow = c(1,1))
  # Plot heatmap of count matrix
  select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
  df <- as.data.frame(colData(dds)[,c("condition","allele")])
  pheatmap::pheatmap(assay(rld)[select,], cluster_rows=FALSE,cluster_cols=FALSE, annotation_col=df)
  pheatmap::pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
                     cluster_cols=FALSE, annotation_col=df)
  
  # Heatmap of sample to sample dist
  sampleDists <- dist(t(assay(rld)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
  pheatmap::pheatmap(sampleDistMatrix,
                     clustering_distance_rows=sampleDists,
                     clustering_distance_cols=sampleDists,
                     col=colors)
  
  # PCA plot of all samples
  data <- plotPCA(rld, intgroup=c("condition", "allele"), returnData=TRUE)
  percentVar <- round(100 * attr(data, "percentVar"))
  ggplot2::ggplot(data, ggplot2::aes(PC1, PC2, color=condition, shape=allele)) + ggplot2::geom_point(size=3) +
    ggplot2::xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ggplot2::ylab(paste0("PC2: ",percentVar[2],"% variance"))
  
  # Number of diffexp genes
  DESeq2::plotMA(ddr, main="MAplot: Genes with allelic bias", ylim=c(-2,2))
  if(barplot == TRUE){
    annot <- rnaResultObject$annotation[1:2] # merge with ddr to remove chr
    merged  <- merge(as.data.frame(ddr),annot,by.x=0,by.y=1)
    altup <- nrow(dplyr::filter(merged,!(excludeChr %in% Chr), padj < fdrCutoff,log2FoldChange > 1))
    merged <- merge(as.data.frame(ddr),annot,by.x=0,by.y=1)
    refup <-  nrow(dplyr::filter(merged,!(excludeChr %in% Chr),padj < fdrCutoff,log2FoldChange < -1))
    # plot
    rnaResultObject$alleleinfo$refAllele %>% as.character() -> ref
    rnaResultObject$alleleinfo$altAllele %>% as.character() -> alt
    reshape2::melt(data.frame(allele = c(ref,alt), genes = c(refup,altup))) -> plotdat
      ggplot2::ggplot(plotdat,ggplot2::aes(allele,value)) + 
          ggplot2::geom_bar(stat = "identity",position = "dodge") +
              ggplot2::labs(y = "Biased genes") + ggplot2::theme_gray(base_size = 16)
  }
  dev.off()
}

#' Write the output
#'
#' @param rnaResultObject output from alleleDiff_rna function
#' @param annotateFrom either "dataset","ensembl" or NULL. telling where to annotate the output from
#' @param species Species for annotation (if annotateFrom = "ensembl")
#' @param excludeChr any chromosome to exclude from output
#' @param fdrCutoff fdr cutoff (same as in previous script)
#' @param outfileName File to write the output in
#' @return csv file with output
#' @examples
#' writeOutput_rna(rnaResultObject,annotateFrom = "dataset", species = "Mus musculus",
#' excludeChr = "chr12",fdrCutoff = 0.01,outfileName)
#' 

writeOutput_rna <- function(rnaResultObject,annotateFrom = "dataset", species = "Mus musculus",
                            excludeChr = "chr12",fdrCutoff = 0.01,outfileName){
  
  results <- as.data.frame(rnaResultObject$DEresult)
  ## Get annotation
  if(annotateFrom == "ensembl"){
    if(is.null(species)) stop("Need species name to annotate from ensembl")
    message("Fetching annotations from ENSEMBL")
    # Download data from ensembl
    mart <- biomaRt::useMart("ensembl", path="/biomart/martservice")
    dataset <- biomaRt::listDatasets(mart)
    dataset <- as.character(dataset[grep(species,dataset[,2]),1])
    mart <- biomaRt::useDataset(dataset, mart=mart)
    info <- biomaRt::getBM(attributes=c("ensembl_gene_id", "chromosome_name", "strand", "start_position", 
                                        "end_position", "gene_biotype","external_gene_name"), mart=mart)
    info$strand <- gsub("-1","-",info$strand)
    info$strand <- gsub("1","+",info$strand)
  } else if(annotateFrom == "dataset"){
    if(is.null(rnaResultObject$annotation)) stop("No annotation available in dataset")
    message("Annotating from within datset")
    info <- rnaResultObject$annotation
    }else {
      warning("No annotation found. You wont get annotated output")
    }
  # Add annotation
  results <- merge(results,info,by.x = 0,by.y=1,all.x = TRUE)
  
  ## Filter by Cutoff or chromosome if given
  if(!(is.null(fdrCutoff))){
    results <- results[results$padj < fdrCutoff,]
  }
  if(!(is.null(excludeChr))){
    results <- dplyr::filter(results, !(excludeChr %in% Chr))
  }
  
  ## Add another column for allelic bias
  refAllele = rnaResultObject$alleleinfo$refAllele
  altAllele = rnaResultObject$alleleinfo$altAllele
  results$Status <- ifelse(results$log2FoldChange > 1, paste0(altAllele,"_biased"),
                           ifelse(results$log2FoldChange < -1, paste0(refAllele,"_biased"),
                                  "NA"))
  
  ## Write back
  write.csv(results,file = outfileName,quote = FALSE,row.names = FALSE)
  
}
