
#loading packages 
library(org.Mm.eg.db)
library("DESeq2")
library("genefilter")
library(ggplot2)
library(ggrepel)

#Loading Dataset

data <- load('/ihome/yufeihuang/zhl169/timothy/Gao_mouse_CASTOR1/feature_counts.RData')
counts_KO <- count_matrix_KO$counts
counts_WT <- count_matrix_WT$counts
count_matrix <- cbind(counts_KO, counts_WT)

genes <- rownames(count_matrix)
genes <- sapply(strsplit(genes, ".",fixed = TRUE), function(x) x[1])
mapped_id <- mapIds(org.Mm.eg.db,keys = genes,column = 'SYMBOL',keytype = 'ENSEMBL')
na_index <- is.na(mapped_id)
rownames(count_matrix) <- mapped_id
count_matrix <- count_matrix[!na_index,]
colnames(count_matrix) <- sapply(strsplit(colnames(count_matrix), ".",fixed = TRUE), function(x) x[1])

#this will return final data matrix where each row is a gene expression and column is a sample
#Starting from count matrices
count_matrix <- as.data.frame(count_matrix)

#summary of conditions 
cond <- sapply(strsplit(colnames(count_matrix), "_",fixed = TRUE), function(x) x[1])
cond <- sapply(strsplit(cond, "-",fixed = TRUE), function(x) x[2:4])
cond[is.na(cond)] <- '0h'
treatment <- cond[1,]
time <- cond[3,]
cond <- paste0(cond[1,],'_',cond[3,])
colData <- data.frame(row.names=colnames(count_matrix), cond,treatment,time)

#DEseq dataset object from the matrix of counts and the sample information table
dds <- DESeqDataSetFromMatrix(count_matrix, colData, ~ treatment)

#removing rows of the DESeqDataSet that have no counts, or only a single count across all samples.

#32875 genes before filtering
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds) 
#24067 genes after filtering

#Testing suitable transformation choice to stabilize variance and mean. 
#rlog
rld <- rlog(dds, blind = FALSE)

#vst
vsd <- vst(dds, blind = FALSE)

library("dplyr")
library("ggplot2")

dds <- estimateSizeFactors(dds)

colnames(dds) <- c(colData$cond)
colnames(rld) <- c(colData$cond)
colnames(vsd) <- c(colData$cond)

#select columns of samples we want to plot
samples_to_plot <- c('KO1_0h', 'WT1_0h')

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, samples_to_plot]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(rld)[, samples_to_plot]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vsd)[, samples_to_plot]) %>% mutate(transformation = "vst"))

ggplot(df, aes(x = KO1_0h, y = WT1_0h)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation) 

#Sample eucledian distances
sampleDists <- dist(t(assay(rld)))
library("pheatmap")
library("RColorBrewer")

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- c(colnames(rld))
colnames(sampleDistMatrix) <- c(colnames(rld))
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
write.csv(sampleDistMatrix, "/ix/yufeihuang/hasib/CASTOR_GAO/results/samples_eld_distance.csv")

#PCA plot

pcaData <-plotPCA(rld, intgroup = c('treatment', 'time'))
percentVar <- round(100 * attr(pcaData$data, "percentVar"))

ggplot(pcaData$data, aes(x = PC1, y = PC2, color = treatment, shape = time)) +
  geom_point(size =3) +
  xlab("PC1") +
  ylab("PC2") +
  coord_fixed()

#multidimensional scaling (MDS) plot

mds <- as.data.frame(colData(rld))  %>%
  cbind(cmdscale(sampleDistMatrix))

mds <- as.data.frame(colData(rld))
cmd_scale <- cbind(cmdscale(sampleDistMatrix))
colnames(cmd_scale) <- c('x', 'y')
mds_2 <- cbind(mds, cmd_scale)

ggplot(mds_2, aes(x = x, y = y, color = treatment, shape = time)) +
  geom_point(size = 3) + coord_fixed() 

#Suggestion to select samples 
sampleDistMatrix <- as.data.frame(sampleDistMatrix)
temp_Data <- sampleDistMatrix[grepl("_12hr", names( sampleDistMatrix )) , grepl("_12hr", names( sampleDistMatrix ))]

#differential analysis on every time point

result<- data.frame()
for (hasib in unique(colData$time)) {
  temp_Data <- count_matrix
  colnames(temp_Data) <- c(colData$cond)
  temp_Data <- temp_Data[ , grepl( paste("_", hasib, sep = "") , names( temp_Data )) ]
  temp_cond <- colData[grepl( paste("_", hasib, sep = "") , colData$cond), ]
  rownames(temp_cond) <- c(temp_cond$cond)
  temp_cond$treatment <- substr(temp_cond$treatment,1,nchar(temp_cond$treatment)-1)
  dds <- DESeqDataSetFromMatrix(temp_Data, temp_cond , ~ treatment)
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)
  res <- results(dds)
  res2 <- na.omit(res)
  res2 <- res2[res2$padj < 0.05,]
  res_df <- data.frame(res2)
  res_df$condition <- hasib
  res_df$DEG <- rownames(res_df)
  result <- rbind(result , res_df)
  #printing highly variable genes
  vst <- vst(dds, blind = FALSE)
  topVarGenes <- head(order(rowVars(assay(vst)), decreasing = TRUE), 30)
  mat  <- assay(vst)[ topVarGenes, ]
  mat  <- mat - rowMeans(mat)
  anno <- as.data.frame(colData(vst)[, c("cond","treatment")])
  
  pdf(file = paste("/ix/yufeihuang/hasib/CASTOR_GAO/figures/all_samples_", hasib,'.pdf' ,sep = ""))
  pheatmap(mat, annotation_col = anno)
  dev.off()
}
write.csv(result, "/ix/yufeihuang/hasib/CASTOR_GAO/results/DEG_all_samples.csv")

ggplot(data=result, aes(x=log2FoldChange, y=-log10(padj), col=DEG, label=condition)) + 
  geom_point(size = 0.3) + theme_minimal() + geom_text_repel(aes(label = DEG), size = 3, show.legend = FALSE)