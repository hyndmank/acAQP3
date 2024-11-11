#### Load libraries
library(DESeq2)
library(ggplot2)
library(readr)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(gplots)
library(calibrate)
library(AnnotationDbi)
library(EnsDb.Mmusculus.v79)
library (vsn)
library(EnhancedVolcano)
library(biomaRt)
library(dplyr)
library(stringr)

############################ SETTING UP ################################
#### Load data

## Import data
countdata <- as.matrix(read.table("counts_clean.txt", header = TRUE, row.names=1))
metadata <- read.csv("metadata_new1.csv") # created 6 groups (MWT, MQ, MR, FWT, FQ, FR)

## Remove .bam or .sam from column names
colnames(countdata)<-gsub("\\Aligned.sortedByCoord.out.bam$", "", colnames(countdata))

#### Set up BiomaRT for gene annotation 

## Set up ensembl parameters
ensembl <- useEnsembl(biomart = 'genes', 
                      dataset = 'mmusculus_gene_ensembl', # choose dataset you want from the available data sets
                      version = 103)

## Set the filter type and values
ourFilterType <- "ensembl_gene_id"

## Set the list of attributes (choose any of the attributes you like from the list ensembl_attributes)
attributeNames <- c("ensembl_gene_id",
                    "external_gene_name", 
                    "entrezgene_id", 
                    "description",
                   "chromosome_name", 
                    "start_position",
                    "end_position",
                    "gene_biotype", 
                    "pdb",
                   "phenotype_description")

########################## DESEQ2 ##################################

#### Construct DESEQDataSet Object
dds <- DESeqDataSetFromMatrix(countData=countdata, 
                              colData=metadata, 
                              design=~Group)
## prefilter low count genes
keep <- rowSums(counts(dds)) >= 100
dds <- dds[keep,]

## Run DESeq
dds <- DESeq(dds) 
resultsNames(dds)

#### QUality control
## Log transformation
rld <- rlogTransformation(dds, blind = TRUE)
head(rld)

## Plot dispersion
plotDispEsts(dds, main = "Dispersion plot")

## heatmap and sample to sample distance
sampleDists <- dist( t( assay(rld) ) )
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

## Plot PCA
plotPCA(rld, intgroup = "Group")

## this gives log2(n + 1)
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))
vsd <- vst(dds, blind=FALSE)
meanSdPlot(assay(vsd))

## heatmap of samples
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Group", "ID")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

## look for sample outliers
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)


#### Pairwise comparisons

## Create a list for looping
# Male WT vs Q
M_Q.WT <- c("MQ", "MWT", "Male Q vs WT", "M_Q.WT")
# Male WT vs R
M_R.WT <- c("MR", "MWT", "Male R vs WT", "M_R.WT")
# Male Q vs R
M_Q.R <- c("MQ", "MR", "Male Q vs R", "M_Q.R")
# Female WT vs Q
F_Q.WT <- c("FQ", "FWT", "Female Q vs WT", "F_Q.WT")
# Female WT vs R
F_R.WT <- c("FR", "FWT" ,"Female R vs WT", "F_R.WT")
# Female Q vs R
F_Q.R <- c("FQ", "FR", "Female Q vs R", "F_Q.R")
# WT Female vs Male 
WT_F.M <- c("FWT", "MWT", "WT Female vs Male", "WT_F.M")
# Q Female vs Male
Q_F.M <- c("FQ", "MQ", "Q Female vs Male", "Q_F.M")
# R Female vs Male
R_F.M <- c("FR", "MR", "R Female vs Male", "R_F.M")


my_list <- list(M_Q.WT, M_R.WT, M_Q.R, F_Q.WT, F_R.WT, F_Q.R, WT_F.M, Q_F.M, R_F.M)

#### Loop
lapply(my_list, function(i){
  ## Run results
  res <- results(dds, contrast=c("Group", i[1], i[2]), alpha = 0.05)
  head(res)
  summary(res)
  
    ## Examine plot of p-values
  file_name <- paste(i[3],"_histogram_pvalue.png", sep="")
  png(file_name, 1000, 1000, pointsize = 20)
  hist(res$pvalue, breaks=50, col="grey")
  dev.off()
  
  ## Examine independent filtering
  file_name <- paste(i[3],"_quantiles_of_BaseMean.png", sep="")
  png(file_name, 1000, 1000, pointsize = 20)
  metadata(res)$filterThreshold  # outdated codes: attr(res, "filterThreshold")
  plot(metadata(res)$filterNumRej, type="b", xlab="quantiles of baseMean", ylab="number of rejections") # outdated codes: attr(res,"filterNumRej")
  dev.off()
  
  ## get the Ensembl IDs from our results table
  filterValues <- rownames(res)
  
  ## run the annotation query
  annot <- getBM(attributes=attributeNames, 
                 filters = ourFilterType, 
                 values = filterValues, 
                 mart = ensembl)
  head(annot)
  
  ## Merge with read counts of all samples
  restable <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  
  ## Annotate Ensembl ID using annot table
  names(restable)[1] <- "EnsemblID"
  restable$GeneName <- with(restable, annot$external_gene_name[match(EnsemblID, annot$ensembl_gene_id)])
  restable$EntrezID <- with(restable, annot$entrezgene_id[match(EnsemblID, annot$ensembl_gene_id)])
  restable$Chromosome <- with(restable, annot$chromosome_name[match(EnsemblID, annot$ensembl_gene_id)])
  restable$Start <- with(restable, annot$start_position[match(EnsemblID, annot$ensembl_gene_id)])
  restable$End <- with(restable, annot$end_position[match(EnsemblID, annot$ensembl_gene_id)])
  restable$GeneType <- with(restable, annot$gene_biotype[match(EnsemblID, annot$ensembl_gene_id)])
  restable$Pdb <- with(restable, annot$pdb[match(EnsemblID, annot$ensembl_gene_id)])
  restable$Phenotype <- with(restable, annot$phenotype_description[match(EnsemblID, annot$ensembl_gene_id)])
  
  ## Reorder columns
  restable <- restable %>% relocate(GeneName, .after = EnsemblID)
  restable <- restable %>% relocate(EntrezID, .after = GeneName)
  restable <- restable %>% relocate(Chromosome, .after = EntrezID)
  restable <- restable %>% relocate(Start, .after = Chromosome)
  restable <- restable %>% relocate(End, .after = Start)
  restable <- restable %>% relocate(GeneType, .after = End)
  restable <- restable %>% relocate(Pdb, .after = GeneType)
  restable <- restable %>% relocate(Phenotype, .after = Pdb)
  
  ## Sort the results based on adjusted p-values (smallest to largest) 
  restable <- restable %>% dplyr::arrange(padj)
  restable <- as.data.frame(restable)
  
  ## Save restable object
  filename <- paste(i[4],".rds", sep="")
  saveRDS(restable, filename)
  
  ## plot most variable gene
  plotCounts(dds, gene=which.min(res$padj), intgroup=c("Group"), pch = 19)
  
  summary(restable)
  
  ## Volcano plot with "significant" genes labeled
  volcanoplot <- function (res, lfcthresh=0.5, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1.5, ...) {
    with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
    with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
    with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="blue", ...))
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
    if (labelsig) {
      require(calibrate)
      with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=GeneName, cex=textcx, ...))
    }
    legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("green","blue","red"))
  }
  file_name <- paste(i[3],"_DE_volcanoplot.png", sep="")
  png(file_name, 1200, 1000, pointsize=20)
  volcanoplot(restable, lfcthresh=0.5, sigthresh=0.05, textcx=1.5, xlim=c(-2.3, 2))
  dev.off()

  ## MA plot
  maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
    with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
    with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
    if (labelsig) {
      require(calibrate)
      with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=GeneName, cex=textcx, col=2))
    }
  }
  file_name <- paste(i[3],"_DE_MAplot.png", sep="")
  png(file_name, 1500, 1000, pointsize=20)
  maplot(restable, main="MA Plot")
  dev.off()
  
  ## Write results
  file_name <- paste(i[3],"_DE_Results.csv", sep="")
  write.csv(restable , file=file_name)

})
