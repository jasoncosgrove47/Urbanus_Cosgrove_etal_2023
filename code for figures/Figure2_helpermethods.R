##########################################################
# This file contains a set of helper methods that 
# are used in the BarcodeMouse10X.Rmd script
#
# @author: Jason Cosgrove (jason.cosgrove47@gmail.com)
# @date: 05/02/2020
#
##########################################################

# this method returns a vector of booleans
# with TRUE if an element of x is found in y
'%ni%' <- function(x,y)!('%in%'(x,y))



geneSetAnnotation <- function(lsks, gene.set.names){
  gene.sets <- read.csv("genesets.csv")
  for(i in 1:length(gene.set.names)){
    genes.of.interest <- paste(gene.sets[,gene.set.names[i]])
    genes.of.interest <- genes.of.interest[genes.of.interest != ""]
    genes.of.interest <- intersect(genes.of.interest,rownames(lsks@assays$SCT))
    lsks <- AddModuleScore(lsks, features = list(genes.of.interest),name = gene.set.names[i])
  }
  return(lsks)
  
}


#load in the outputs of cell-ranger, and convert into a seurat object adding some metadata
generateSeuratObject <- function(){
  
  #load in the outputs of the cellranger pre-processing pipeline
  gene_bc_matrix.neg <- Read10X(data.dir="GFPneg") 
  gene_bc_matrix.pos <-Read10X(data.dir= "GFPpos") 
  
  dataset <- cbind(gene_bc_matrix.pos,gene_bc_matrix.neg)
  colnames(dataset) <- make.unique(colnames(dataset), sep = "_")
  
  #generate a dataframe with the batchID for each cell
  Batch <- c(rep("GFPpos",642),rep("GFPneg",2231))
  cell_anns <- data.frame(batch = Batch)
  rownames(cell_anns) <- colnames(dataset)
  
  print(paste("dataset contains" , nrow(dataset),"genes"), sep = " " )
  print(paste("dataset contains" , ncol(dataset),"cells"), sep = " " )
  
  #now we have the dataset and the metadata we can make the seurat object
  lsks <- CreateSeuratObject(counts = dataset, min.cells = 10, min.features = 500,meta.data = cell_anns, 
                             project = "barcode")
  
  return(lsks)
}


# this method converts a seurat object into a single cell experiment object
# and produces PCA plots overlaying information about key quality control metrics
# like library size
checkNormalisation <- function(assay,lsks){
  if(assay == "SCT"){
    lsks <- RunPCA(lsks, npcs = 50, verbose = FALSE,assay = "SCT")
    lsks.sce <- as.SingleCellExperiment(lsks,assay = "SCT")
  }
  else if(assay == "RNA"){
    lsks <- RunPCA(lsks, npcs = 50, verbose = FALSE,assay = "RNA")
    lsks.sce <- as.SingleCellExperiment(lsks,assay = "RNA")
  }
  
  p1 <- plotPCASCE(
    lsks.sce,
    colour_by = "nCount_RNA",
    size_by="nFeature_RNA",
    by_exprs_values= "counts"
    #rerun=TRUE,
    #run_args = list(by_exprs_values= "counts")
  ) + ggtitle("raw data")
  
  p2 <- plotPCASCE(
    lsks.sce,
    colour_by = paste("nCount_",assay, sep = ""),
    size_by=paste("nFeature_",assay, sep = ""),
    by_exprs_values= "logcounts"
    #rerun=TRUE,
    #run_args = list(by_exprs_values= "logcounts")
  ) + ggtitle(paste("normalised_",assay, sep = ""))
  
  return(list(p1 = p1, p2 = p2))
  
}
