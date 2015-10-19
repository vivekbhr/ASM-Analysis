
# Functions to perform allele-sp differential expression
# (c) Vivek Bhardwaj (bhardwaj@ie-freiburg.mpg.de)
# First created: 13th May 2015

### Read the Files and Count features for the RNA-Seq Samples

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
    fcres <- Rsubread::featureCounts(files = bam.files,annot.ext = AnnotationFile,
                                     isGTFAnnotationFile=TRUE,...)
  } else {
    message("Counting reads using standard parameters")
    fcres <- Rsubread::featureCounts(files = bam.files,annot.ext = AnnotationFile,
                                     isGTFAnnotationFile=TRUE,useMetaFeatures=TRUE,
                                     isPairedEnd=TRUE,requireBothEndsMapped=FALSE,checkFragLength=FALSE,
                                     nthreads=nthreads,strandSpecific=2,minMQS=0,countPrimaryAlignmentsOnly=TRUE)
  }
  # Save Fcount output 
  colnames(fcres$counts) <- samp[,2] # change colnames to samplenames
  fcres$targets <- samp[,2] # changed targets(just another df) to samplenames
  if(!(is.null(outfileName))){
    write.table(fcres$counts,outfileName,sep="\t",row.names=F,quote=F)
  }
  
  ## Make the output
  fcres$design <- design
  fcres$alleleinfo <- alleleinfo
  return(fcres)
}

### Run DESeq2 with interaction design

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
  return(rnaCountObject)
}

### Write the output

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
    results <- results[results$Chr != excludeChr,]
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
