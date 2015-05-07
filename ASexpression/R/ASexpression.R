
# library "ASexpression"
# Author: Vivek Bhardwaj, MPI-IE
# Date: 7th May 2015


readFCres <- function(CountFile,ContSampleNames,TestSampleNames){
  fc.asgenes <- read.table(CountFile,header=T,sep="\t",row.names= 1)
  contNames <- as.character(read.delim(ContSampleNames,sep="\n",header=F)$V1) # samplenames should be given as a file
  testNames <- as.character(read.delim(TestSampleNames,sep="\n",header=F)$V1) # samplenames should be given as a file
  fcres <- list(cont = fc.asgenes[grep(paste(contNames,collapse = "|"),colnames(fc.asgenes))],
                test = fc.asgenes[grep(paste(testNames,collapse = "|"),colnames(fc.asgenes))])  
  return(fcres)
}

runDESeq <- function(fcres,design.cont,design.test){ #design = c(MatOverPat,PatOverMat)
  library('DESeq2')
  contMatrix <- design.cont
  testMatrix <- design.test
  
  contMatrix$design = paste(design.cont$sample,design.cont$allele,sep="_")
  testMatrix$design = paste(design.test$sample,design.test$allele,sep="_")
  
  contMatrix$design <- as.factor(contMatrix$design)
  testMatrix$design <- as.factor(testMatrix$design)
  sampList <- list(cont = DESeqDataSetFromMatrix(countData = as.matrix(fcres$cont),colData = contMatrix, design = ~ design),
                   test = DESeqDataSetFromMatrix(countData = as.matrix(fcres$test),colData = testMatrix, design = ~ design))
  sampList <- list(cont = DESeq(sampList$cont), test= DESeq(sampList$test))
  return(sampList)
}

pullResults <- function(DESeqOutputList,maternal="129S1",paternal="CASTEiJ",
                        compareGroup = "MatOverPat", design.cont, design.test){
  design.cont$design = paste(design.cont$sample,design.cont$allele,sep="_")
  design.test$design = paste(design.test$sample,design.test$allele,sep="_")
  
  biglist <- list(cont = design.cont[!duplicated(design.cont$design),],
                  test = design.test[!duplicated(design.test$design),])
  
  output <- list()
  for(group in names(biglist)){
    output[[group]] <- list()
    for(name in unique(biglist[[group]]$sample)){
      matSample = paste(name,maternal,sep="_")
      patSample = paste(name,paternal,sep="_")
      if(compareGroup == "MatOverPat") testSet <- c(matSample,patSample) else testSet <- c(patSample,matSample)
      contrasts =  c("design",testSet)
      output[[group]][[name]] = results(DESeqOutputList[[group]],contrast = contrasts, cooksCutoff = FALSE,alpha = 0.01)
    }
  }
  return(output)
}

filterByControls <- function(DESeqResultList,matchfile,padjCutoff = 0.01){
  diff.unique = list()
  for(tf in names(DESeqResultList$test)){
    matchCont <- as.character(matchfile[which(matchfile$test == tf),'controls'])
    diffgenes <- subset(as.data.frame(DESeqResultList$test[[tf]]), padj < padjCutoff)
    diffgenes.cont <- subset(as.data.frame(DESeqResultList$cont[[matchCont]]),padj < padjCutoff)
    uniqueToTF <- setdiff(row.names(diffgenes),row.names(diffgenes.cont))
    diff.unique[[tf]] <- diffgenes[which(row.names(diffgenes) %in% uniqueToTF),]
  }
  return(diff.unique)
}

makeSomePlots <- function(DESeqOutputList,ControlFilteredOutput,compareGroup = "MatOverPat"){
  #plot PCA
  library(gridExtra)
  library(plyr)
  library(reshape)
  library(ggplot2)
  PCA <- list()
  for(n in names(DESeqOutputList)){
    PCA[[n]] <- plotPCA(DESeqOutputList[[n]],intgroup = "design")
  }
  grid.arrange(PCA[[1]],PCA[[2]],ncol=2, main= "PCA of samples")
  
  # plot The filtered output stats
  stats <- list()
  for(tf in names(ControlFilteredOutput)){
    if(compareGroup == "MatOverPat"){
      stats[[tf]]$Maternal_Biased <- length(which(ControlFilteredOutput[[tf]]$log2FoldChange > 0)) 
      stats[[tf]]$Paternal_Biased <- length(which(ControlFilteredOutput[[tf]]$log2FoldChange < 0))  
    }else{
      stats[[tf]]$Maternal_Biased <- length(which(ControlFilteredOutput[[tf]]$log2FoldChange < 0)) 
      stats[[tf]]$Paternal_Biased <- length(which(ControlFilteredOutput[[tf]]$log2FoldChange > 0))
    }
  }
  stats <- ldply(stats,as.data.frame)
  stats <- melt(stats)
  ggplot(stats,aes(.id,value,fill=variable,group=variable)) + geom_bar(stat = "identity",position = "dodge") +
    labs(x = "Sample", y = "No. of Genes",fill = "Status",title = "No. of Genes with Change in Allelic Expression")
}

pathwayEnrichment <- function(ControlFilteredOutput,organism="mmu"){
  GOSeqInput = list()
  for(n in names(ControlFilteredOutput)){
    undiff <-rep(0,length(setdiff(rownames(ase.results$test[[n]]),rownames(ControlFilteredOutput[[n]])) ))
    names(undiff) <- setdiff(rownames(ase.results$test[[n]]),rownames(ControlFilteredOutput[[n]]))
    # had to create another vector undiff and add to the input, as GOSeq needs all genes as input  
    GOSeqInput[[n]]$up <- as.integer(ControlFilteredOutput[[n]]$log2FoldChange > 0)
    names(GOSeqInput[[n]]$up) <- rownames(ControlFilteredOutput[[n]])
    GOSeqInput[[n]]$up <- c(GOSeqInput[[n]]$up,undiff)
    GOSeqInput[[n]]$down <- as.integer(ControlFilteredOutput[[n]]$log2FoldChange < 0)
    names(GOSeqInput[[n]]$down) <- rownames(ControlFilteredOutput[[n]])
    GOSeqInput[[n]]$down <- c(GOSeqInput[[n]]$down,undiff)
  }
  
  
  # Using GOSeq for pathway enrichment analysis
  library('goseq')
  library('KEGGREST')
  
  keggid2name <- keggList("pathway", "mmu")
  names(keggid2name) <- sapply(names(keggid2name), substring, 9)
  
  GOSEQ <- list()
  for(n in names(GOSeqInput)){
    for(cat in c('up','down')){
      null_model <- nullp(GOSeqInput[[n]][[cat]], genome= "mm9",id="ensGene")
      kegg.enrich <- goseq(null_model,genome = "mm9",id="ensGene",test.cats= c("KEGG"),method = "Wallenius")
      kegg.fdrs <- p.adjust(kegg.enrich[, 2], method = "BH")
      kegg.enrich <- cbind(kegg.enrich, FDR = kegg.fdrs)
      kegg.enrich <- cbind(kegg.enrich, pathway = keggid2name[kegg.enrich$category])
      GOSEQ[[n]][[cat]] <- kegg.enrich[kegg.enrich$FDR < 0.05,]    
    }
  }
  
  GOSEQ <- lapply(GOSEQ, function(x) x = ldply(x,data.frame))
  return(GOSEQ)
}
