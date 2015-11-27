
#### ~~~~ Functions to Run CSAW as part of AS analysis pipeline ~~~~ ####
### (c) Vivek Bhardwaj (bhardwaj@ie-freiburg.mpg.de)

#' Read the Files and Count windows for ChIP-Seq Samples
#'
#' @param csvFile csvfile with sample information
#' @param refAllele Reference Allele (name must match the name in samplesheet)
#' @return chipCountObject : a list with window counts and sampleinfo
#' @examples
#' readfiles_chip(csvFile = "testBAMs/testSampleSheet.csv", refAllele = "pat")
#' 
readfiles_chip <- function(csvFile = "testBAMs/testSampleSheet.csv", refAllele = "pat"){
  
  # Parse Sample sheet
  samp <- read.csv(csvFile,header=TRUE, 
                   colClasses = c("character","character","factor","factor","factor","character"))
  if(any(grepl("chip",samp[,1]))){
    message("Extracting sample Info")
    samp <- samp[which(grepl("chip",samp[,1])),]
  } else stop("No Chip Samples found for diffBinding test. Check sample sheet")
  
  # readFiles using CSAW
  bam.files <- samp[,6]
  pe.param <- csaw::readParam(max.frag=500, pe="both") # use the param for pe reads
  message("Counting reads in windows")
  counts <- csaw::windowCounts(bam.files = bam.files,param=pe.param)
  
  # make model matrix for edgeR
  design <- data.frame(row.names = samp[,2] , samp[,3:5]) ## Add: A warning if rownames not uniq
  colnames(design) <- c("condition","allele","tf")
  design[,2] <- relevel(design[,2],refAllele)
  design[,3] <- droplevels(design[,3])
  design <- model.matrix(~ allele * tf,data = design) ## bug : uses colnames to create model matrix!
  # output
  chipCountObject <- list(windowCounts = counts, design = design, sampledata = samp)
  return(chipCountObject)
}

#' Make plots to select window size and pe-distance cutoffs
#'
#' @param chipCountObject output from readfiles_chip
#' @param outdir output directory to write the plots
#' @return chipCountObject : a list with window counts and sampleinfo
#' @examples
#' makeQCplots_chip(chipCountObject,outdir)
#' 

makeQCplots_chip <- function(chipCountObject,outdir){
  
  ## Parse Sample data
  samp <- chipCountObject$sampledata
  bam.files <- samp[,6]
  
  ## Histogram to check frag size cutoff
  message("Checking fragment sizes")
  pdf(paste0(outdir,"/pefrag-sizes.pdf"))
  out <- list(sapply(bam.files,function(x) csaw::getPESizes(x)))
  #plot (change the number of plots)
  for(i in seq(1,16,2)) hist(out[[1]][[i]], breaks=50, xlab="Fragment sizes (bp)",
                             ylab="Frequency", main=paste0("sample_",i),col="steelblue") %>% 
    abline(v=400, col="red",lwd = 3)
  dev.off()
  
  ## Checking cross-correlation
  message("Checking cross-correlation")
  max.delay <- 500
  dedup.on <- csaw::readParam(dedup=TRUE, minq=20)
  x <- csaw::correlateReads(bam.files, max.delay, param=dedup.on)
  #plot
  pdf(paste0(outdir,"/peCross-correlation.pdf"))
  plot(0:max.delay, x, type="l", ylab="CCF", xlab="Delay (bp)")
  dev.off()
  
  ## Choosing appropriate window size
  message("Checking appropriate window sizes")
  plotwc <- function(curbam){
    #print(curbam)
    pe.param <- csaw::readParam(max.frag=500, pe="both") # use the param for pe reads
    windowed <- csaw::windowCounts(curbam, spacing=50, width=50, param= pe.param, filter=20)
    rwsms <- rowSums(assay(windowed))
    maxed <- csaw::findMaxima(rowRanges(windowed), range=1000, metric=rwsms)
    curbam.out <- csaw::profileSites(curbam, rowRanges(windowed)[maxed],
                                        param=pe.param, weight=1/rwsms[maxed])
    return(curbam.out)
  }
  collected <- lapply(bam.files,plotwc) # only plotted 4 bam.files[c(1,3,5,7)] for example
  xranged <- as.integer(names(collected[[1]]))
  # plot
  pdf(paste0(outdir,"/peWindow-sizes.pdf"))
  plot(xranged, collected[[1]], type="l", col="blue", xlim=c(-1000, 1000), lwd=2,
       xlab="Distance (bp)", ylab="Relative coverage per base")
  lines(xranged, collected[[2]], col="forestgreen", lwd=2)
  lines(xranged, collected[[3]], col="grey20", lwd=2)
  lines(xranged, collected[[4]], col="purple", lwd=2)
  legend("topright", col=c("blue", "forestgreen","grey20","purple"),
         paste0("Sample_",1:4), pch=16)
  abline(v=c(-150,200), col="dodgerblue", lty=2)
  dev.off()
  
}

#' Filtering using Input windows
#'
#' @param chipCountObject output from readfiles_chip
#' @param priorCount Minimum count cutoff for windows
#' @return Filtered chipCountObject
#' @examples
#' filterByInput_chip(chipCountObject,priorCount = 5)
#' 

filterByInput_chip <- function(chipCountObject,priorCount = 5){
        
        # Parse Sample data
        samp <- chipCountObject$sampledata
        countdat <- chipCountObject$windowCounts
        # Seperate control and chip samples
        control <- samp[which(samp[,3] == "control"),6]
        control <- countdat[,which(colData(countdat)$bam.files %in% control)]
        chip <- samp[which(samp[,3] == "test"),6]
        chip <- countdat[,which(colData(countdat)$bam.files %in% chip)]
        # Filter chip by control counts
        filter.stat <- csaw::filterWindows(chip, control, type="control", prior.count = priorCount) # min count in window should be 5
        keep <- filter.stat$filter > log2(3) 
        countdat <- countdat[keep,] # now "countdat" contains both input and chip counts, but filtered
        
        # return the same chipCountObject back, but this time filtered
        chipCountObject$windowCounts <- countdat
        return(chipCountObject)
}


### 
#' TMM normalize (get the normfactors out) using given window size
#'
#' @param chipCountObject output from filterByInput_chip
#' @param binsize Size of bins to calculate the normalization factors
#' @param plotfile file with output plots
#' @return Normalized chipCountObject
#' @examples
#' tmmNormalize_chip(chipCountObject,binsize = 10000, plotfile = "TMM_normalizedCounts.pdf")
#' 

tmmNormalize_chip <- function(chipCountObject,binsize = 10000, plotfile = "TMM_normalizedCounts.pdf"){
        
        # Parse Sample data
        samp <- chipCountObject$sampledata
        bam.files <- samp[,6]
        # Get norm factors
        demo <- csaw::windowCounts(bam.files, bin=TRUE, width = binsize)
        normfacs <- csaw::normOffsets(demo)
        chipCountObject$normFactors <- normfacs
        
        # plot normalized counts
        pdf(plotfile)
        par(mfrow=c(1, 3), mar=c(5, 4, 2, 1.5))
        adj.counts <- edgeR::cpm(csaw::asDGEList(demo), log=TRUE)
        for (i in 1:(length(bam.files)-1)) {
                cur.x <- adj.counts[,1]
                cur.y <- adj.counts[,1+i]
                smoothScatter(x=(cur.x+cur.y)/2+6*log2(10), y=cur.x-cur.y,xlab="A", 
                              ylab="M", main=paste("1 vs", i+1))
                all.dist <- diff(log2(normfacs[c(i+1, 1)]))
                abline(h=all.dist, col="red")
        }
        ## MDS plot to check for replicate variability
        for (top in c(100, 500, 1000, 5000)) {
                limma::plotMDS(adj.counts, main=top, col= as.numeric(samp[,4]),labels=samp[,2], top=top)
        }
        dev.off()
        
        ## Return normfactors
        return(chipCountObject)
}


### 
#' Test for Diff Bound windows using EdgeR (then merge windows into regions)
#'
#' @param chipCountObject output from tmmNormalize_chip
#' @param plotfile file with output plots
#' @param tfname which TF to extract results for (must match with the name in samplesheet)
#' @return chipResultObject with differentially bound regions
#' @examples
#' getDBregions_chip(chipCountObject,plotfile = NULL, tfname = "msl2")
#' 

getDBregions_chip <- function(chipCountObject,plotfile = NULL, tfname = "msl2"){
        
        # Make DGElist
        y <- csaw::asDGEList(chipCountObject$windowCounts, norm.factors = chipCountObject$normFactors)
        design <- chipCountObject$design
        # Estimate dispersions
        y <- edgeR::estimateDisp(y, design)
        o <- order(y$AveLogCPM)
        fit <- edgeR::glmQLFit(y, design, robust=TRUE)
        # and plot dispersions
        if(!(is.null(plotfile))){
                pdf(plotfile)
                par(mfrow=c(1,2))
                plot(y$AveLogCPM[o], sqrt(y$trended.dispersion[o]), type="l", lwd=2,
                     ylim=c(0, 1), xlab=expression("Ave."~Log[2]~"CPM"),
                     ylab=("Biological coefficient of variation"))
                edgeR::plotQLDisp(fit)
                dev.off()
        }
        
        # TEST for DB windows
        coef_one <- colnames(fit)[2] # it's the first coefficient that can be extracted (this will be the non-reference allele)
        tf <- colnames(chipCountObject$sampledata)[5]
        ## other coefs will be added above this, the names of 2nd coef have to be provided as tfname
        results <- edgeR::glmQLFTest(fit, coef = paste0(coef_one,":",tf,tfname))
        
        # Clustering DB windows into regions: Using quick and dirty method
        merged <- csaw::mergeWindows(rowRanges(chipCountObject$windowCounts), tol=500L)
        tabcom <- csaw::combineTests(merged$id, results$table) # get combined test p-value for merged windows
        # Return all results
        chipResultObject <- list(fit = fit, results = results, mergedRegions = merged, combinedPvalues = tabcom)
        return(chipResultObject)
}


### 
#' Annotate and print the output regions
#'
#' @param chipResultObject output from getDBregions_chip
#' @param outfileName name of output files
#' @param annotation TRUE if you want to annotate the regions
#' @param Txdb Txdb object to annotate the output (if annotation = TRUE)
#' @param Orgdb orgdb object to annotate the output (if annotation = TRUE)
#' @param tfname which TF to extract results for (must match with the name in samplesheet)
#' @return File with differentially bound regions
#' @examples
#' writeOutput_chip(chipResultObject, outfileName, annotation = TRUE, 
#' Txdb = TxDb.Mmusculus.UCSC.mm9.knownGene, Orgdb = org.Mm.eg.db)
#' 

writeOutput_chip <- function(chipResultObject, outfileName, annotation = TRUE, fdr = 0.05,
                             Txdb = TxDb.Mmusculus.UCSC.mm9.knownGene, Orgdb = org.Mm.eg.db){
        # get merged regions
        merged <- chipResultObject$mergedRegions
        tabcom <- chipResultObject$combinedPvalues
        ## adding gene annotation
        if(annotation){
                anno <- csaw::detailRanges(merged$region, txdb= Txdb,
                                     orgdb=org.Mm.eg.db, promoter=c(1000, 1000), dist=5000)
                anno.ranges <- csaw::detailRanges(txdb=Txdb, orgdb = Orgdb)
                
                ## Print regions and genes as output
                ofile <- gzfile(paste0(outfileName,".gz"), open="w")
                write.table(as.data.frame(merged$region)[,1:3], tabcom, anno,
                            file=ofile, row.names=FALSE, quote=FALSE, sep="\t")
                close(ofile)
                
        } else {
                ## Print regions as bed files
                is.sig <- tabcom$FDR <= fdr
                test <- merged$region[is.sig]
                if(length(test) > 0){
                  test$score <- -10*log10(tabcom$FDR[is.sig])
                  names(test) <- paste0("region", 1:sum(is.sig))
                  export(test, paste0(outfileName,".bed"))
                } else {
                  message("output empty! please lower the fdr threshold.")
                }
                
                
        }
        
}

