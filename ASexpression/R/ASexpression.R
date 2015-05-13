
# library "ASexpression"
# Author: Vivek Bhardwaj, MPI-IE
# Date: 13th May 2015


readFCres <- function(CountFile,SampleNames){
  fc.asgenes <- read.table(CountFile,header=T,sep="\t",row.names= 1)
  
  print("Separating samples by Control-Test pairs.")
  fcres <- list()
  for(i in 1:nrow(SampleNames)){
    tf = as.character(SampleNames[i,1])
    set <- as.character(t(SampleNames[i,]))
    fcres[[tf]] <- fc.asgenes[grep(paste(set,collapse = "|"),colnames(fc.asgenes))]
  }
  return(fcres)
}

createDesignMatrix <- function(fcres, SampleNames, baseAllele = "CASTEiJ"){
  design <- list()
  i = 0
  for(n in names(fcres)){
    i = i+1
    rownames <- colnames(fcres[[n]])
    condition <- gsub('(.*)_[0-9]_(.*)','\\1',colnames(fcres[[n]]))
    allele <- as.factor(gsub('(.*)_[0-9]_(.*)','\\2',colnames(fcres[[n]])))
    Control = as.character(SampleNames[i,2])
    
    design[[n]] <- data.frame(row.names = rownames, 
                              condition = relevel(as.factor(condition),Control), 
                              allele = relevel(as.factor(allele),baseAllele))
  }
  return(design)
}

runDESeq <- function(fcres,design, baseAllele = "CASTEiJ",topAllele = "129S1",fdrCutoff=0.01,
                     autodesigned = FALSE, SampleNames = NULL){ 
  source("http://bioconductor.org/biocLite.R")
  if(!(require('DESeq2'))) biocLite('DESeq2')
  library('DESeq2')
  if(autodesigned == FALSE){
    if(is.data.frame(design)){
      listname = as.character(SampleNames[1,1])
      design = list(name = design)
      names(design) = listname
    } 
    print("Using external matrix")
    print("setting factors and levels for DESeq")
    i = 0
    for(n in names(design)){
      i = i +1
      colnames(design[[n]]) = c("condition","allele")
      Control = as.character(SampleNames[i,2])
      design[[n]]$condition = relevel(as.factor(design[[n]]$condition),Control) 
      design[[n]]$allele = relevel(as.factor(design[[n]]$allele),baseAllele)
    } 
  } else if(autodesigned == TRUE) {
    print("Using autodesigned matrix")
  }else stop("Please choose TRUE or FALSE")
  
  
  sampList = list()
  for(name in names(fcres)){
    sampList[[name]] <- DESeqDataSetFromMatrix(fcres[[name]],design[[name]],design = ~ condition + allele + condition:allele)
    sampList[[name]] <- DESeq(sampList[[name]],betaPrior = FALSE)  
  }
  
  print("DESeq finished. Now extracting results.")
  output <- list()
  for(name in names(sampList)){
    query = paste0("condition",name,".allele",topAllele)
    output[[name]] = results(sampList[[name]], name = query, cooksCutoff = FALSE,alpha = fdrCutoff)
  }
  
  allOut <- list(dataSet = sampList, Results = output)
  return(allOut)
}

makeSomePlots <- function(DESeqOutputList,baseAllele = "CASTEiJ",topAllele = "129S1",fdrCutoff=0.01){
  #plot PCA
  #if(!(require('gridExtra'))) install.packages('gridExtra')
  if(!(require('plyr'))) install.packages('plyr')
  if(!(require('reshape'))) install.packages('reshape')
  if(!(require('ggplot2'))) install.packages('ggplot2')
  
  #library(gridExtra)
  library(plyr)
  library(reshape)
  library(ggplot2)
  
  PCA <- list()
  for(n in names(DESeqOutputList$dataSet)){
    PCA[[n]] <- plotPCA(DESeqOutputList$dataSet[[n]],intgroup = c("condition","allele"))
  }
  
  # plot The filtered output stats
  stats <- list()
  for(tf in names(DESeqOutputList$Results)){
  stats[[tf]] = data.frame(Allele = c(topAllele,baseAllele),
                           Biased.genes = c(nrow(subset(DESeqOutputList$Results[[tf]], log2FoldChange > 0 & padj < 0.01)),
                                            nrow(subset(DESeqOutputList$Results[[tf]], log2FoldChange < 0 & padj < 0.01)))
  )
  }
  stats <- ldply(stats,as.data.frame)
  plot <- ggplot(stats,aes(.id,Biased.genes,fill=Allele,group=Allele)) + geom_bar(stat = "identity",position = "dodge") +
              labs(x = "Sample", y = "No. of Genes",fill = "Bias towards",title = "No. of Genes with Change in Allelic Expression")
return(plot)
}


pathwayEnrichment <- function(DESeqOutputList,organism="mmu",fdrCutoff=0.01){
  
  print("preparing sample for Input")
  GOSeqInput = list()
  for(n in names(DESeqOutputList$Results)){
    result = DESeqOutputList$Results[[n]] 
    GOSeqInput[[n]]$up <- as.integer(result$padj < fdrCutoff & result$log2FoldChange > 0)
    names(GOSeqInput[[n]]$up) <- rownames(result)
    GOSeqInput[[n]]$up = na.omit(GOSeqInput[[n]]$up)
    GOSeqInput[[n]]$down <- as.integer(result$padj < fdrCutoff & result$log2FoldChange < 0)
    names(GOSeqInput[[n]]$down) <- rownames(result)
    GOSeqInput[[n]]$down = na.omit(GOSeqInput[[n]]$down)
  }
  

  # Using GOSeq for pathway enrichment analysis
  source("http://bioconductor.org/biocLite.R")
  if(!(require('goseq'))) biocLite('goseq')
  if(!(require('KEGGREST'))) biocLite('KEGGREST')
  if(!(require('org.Mm.eg.db'))) biocLite('org.Mm.eg.db')
  
  library('goseq')
  library('KEGGREST')
  print("Running Enrichment Test")
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





