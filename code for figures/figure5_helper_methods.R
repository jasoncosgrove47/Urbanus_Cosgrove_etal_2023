##############################################################
# HELPER METHODS REQUIRED FOR THE ANALYSIS OF VDJ MOUSE
# DATA
#
# @Author: Jason Cosgrove (jason.cosgrove@curie.fr)
# @Date: 02/01/2020
###############################################################


# is X is in the vector Y
# inputs : "X" is an object 
#          "Y" is a vector
# outputs: A boolean stating TRUE if X is in Y
'%ni%' <- function(x,y)!('%in%'(x,y))


preProcessDataset <- function(sobj, nCount_threshold, percent.mt_threshold){
  
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^mt-")
  
  sobj <- subset(sobj, subset = nCount_RNA < nCount_threshold
                 & percent.mt < percent.mt_threshold)
  
  sobj <- removeUninformativeGenes(sobj)
  
  return(sobj)
  
}




removeUninformativeGenes <- function(dataset){
  genes.to.remove <- c(rownames(dataset)[grepl("Gm",rownames(dataset))],
                       rownames(dataset)[grepl("Rik",rownames(dataset))],
                       rownames(dataset)[grepl("Rp",rownames(dataset))])
  
  
  genes  <- rownames(dataset)[rownames(dataset) %ni% genes.to.remove]
  
  dataset <- subset(dataset, features = genes)
  return(dataset)
}





# classify each cell in a scRNAseq dataset into one of 3 cell cycle stages: G1, S, G2M
# using the cyclone package using the ration approach described in Scialdone et al 2015
#
#  inputs:'lsks' A Seurat object
# outputs: A seurat object with cell cycle phases and associated classifier scores stored 
#         in the metadata slot
runCellCycleAnnotation <- function(lsks){
  
  mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
  ensembl <- mapIds(org.Mm.eg.db, keys=rownames(lsks), keytype="SYMBOL", column="ENSEMBL")
  assignments <- cyclone(lsks@assays$integrated@data, mm.pairs, gene.names=ensembl)
  lsks@meta.data$phases <- assignments$phases
  lsks@meta.data$G1_score <- assignments$normalized.scores$G1
  lsks@meta.data$S_score <- assignments$normalized.scores$S
  lsks@meta.data$G2M_score <- assignments$normalized.scores$G2M
  table(lsks@meta.data$phases)
  return(lsks)
}



nearestNeighbourMapping <- function(referenceDataset,queryDataset,genes,plotColor = "purple"){
  # run the PCA
  #pca <- prcomp(t(dahlin@assays$RNA@data[genes,]), tol = 0.1)
  #save(pca,file = "inputs/dahlin_pcrompPCA.Rda")
  #this step is computationally costly so we just load in one we have performed before
  
  load("inputs/dahlin_pcrompPCA.Rda")
  #lets do external positive controls 
  
  #Use the loadings from the ref dataset PCA to get PCA coords for the query dataset
  pca.scaled <- scale(t(queryDataset@assays$RNA@data[genes,]), pca$center, pca$scale) %*% pca$rotation 
  
  #run the nn algo on the PCA coords
  nn <- nn2(pca$x[,1:10],pca.scaled[,1:10], k = 20)
  
  #get the mean umap coords for all nearest neighbours
  umap.coords <- matrix(0, nrow = ncol(queryDataset),ncol = 2)                    
  rownames(umap.coords) <- colnames(queryDataset)
  for(i in 1:ncol(queryDataset)){
    umap.coords[i,] <- colMeans(dahlin@reductions$umap@cell.embeddings[nn$nn.idx[i,],][1:10,])
  }
  
  
  #plot the results
  x <- DimPlot(object = dahlin,reduction = "umap") + NoLegend()
  g <- ggplot_build(x)
  umap.coords <- umap.coords[colnames(queryDataset),]
  
  p <- plot(x$data[,1],x$data[,2],col = "grey 90",pch = 20, cex = 0.5) + 
    points(x = umap.coords[,1], y = umap.coords[,2] ,cex = 1, pch =20,col = plotColor)
  
  return(p)
}


plotCellCycleResults <- function(dataset.integrated){
  df <- (table( dataset.integrated@meta.data$timepoint,dataset.integrated@meta.data$phases))
  
  par(mfrow=c(2,1))
  barplot(t(t(df) /colSums(df)),legend = FALSE, col = c("cornflowerblue", "darkolivegreen3","orange"))
  barplot(df,legend = TRUE, beside = TRUE,col = c("cornflowerblue", "darkolivegreen3","orange"))

}


createDensityPlot <- function(sobj,bins = 100){
  
  tmp.all<-as.data.frame(Embeddings(object = sobj, reduction = "umap"))
  p <- ggplot(tmp.all, aes(x = UMAP_1, y = UMAP_2)) + geom_point(colour="#00000000") + 
    stat_density_2d(aes(fill = stat(level)), geom = "polygon", bins=bins) + 
    scale_fill_gradientn(colors = c("#4169E100","royalblue", "darkolivegreen3","goldenrod1","red")) +
    theme_classic() + 
    theme(legend.position="bottom")
  
  return(p)
}



geneSignatureScoring <- function(sobj,  gene.sets,gene.set.names,assay){
  
  for(i in 1:length(gene.set.names)){
    genes.of.interest <- paste(gene.sets[,gene.set.names[i]])
    genes.of.interest <- genes.of.interest[genes.of.interest != ""]
    genes.of.interest <- intersect(genes.of.interest,rownames(sobj@assays$RNA))
    sobj <- AddModuleScore(sobj , features = list(genes.of.interest),
                           name = gene.set.names[i],replace = TRUE,assay = assay)
  }
  
  return(sobj)
  
}


# This method is uses to visualise the effects of seurats default normalisation approach
# and also the SCT method using PCA based plots
# Inputs: 'assay' the normalisation approach you would like to visualise
# Outputs: A list of ggplot objects
checkNormalisation <- function(assay,seuratobject){
  if(assay == "SCT"){
    lsks.sce <- as.SingleCellExperiment(lsks,assay = "SCT")
  }
  else if(assay == "RNA"){
    lsks.sce <- as.SingleCellExperiment(lsks,assay = "RNA")
  }
  
  p1 <- plotPCASCE(
    lsks.sce,
    colour_by = "nCount_RNA",
    #by_exprs_values= "logcounts",
    size_by="nFeature_RNA",
    rerun=TRUE,
    run_args = list(exprs_values= "counts")
  ) + ggtitle("raw data")
  
  p2 <- plotPCASCE(
    lsks.sce,
    colour_by = paste("nCount_",assay, sep = ""),
    #by_exprs_values= "logcounts",
    size_by=paste("nFeature_",assay, sep = ""),
    rerun=TRUE,
    run_args = list(exprs_values= "logcounts")
  ) + ggtitle(paste("normalised_",assay, sep = ""))
  
  return(list(p1 = p1, p2 = p2))
  
}
