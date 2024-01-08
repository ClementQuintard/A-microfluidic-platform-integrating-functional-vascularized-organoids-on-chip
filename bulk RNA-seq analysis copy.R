# RNA-seq analysis 

# Setup -------------------------------------------------------------------
set.seed(42)
# Packages ----------------------------------------------------------------
library(DESeq2)
# this is only needed for very big analyses
#library("BiocParallel")
#register(MulticoreParam(4))
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(EnhancedVolcano)
library(ggplot2)
library(ggrepel)
library(enrichplot)
library(clusterProfiler)
library(dendextend)
# Data loading ------------------------------------------------------------
# get design matrix
design_matrix <- read.csv("design-matrix-cq.csv", header = T)
# files to extract data from

files <- list.files("Quant genes txt")

# read in first file to grab rownames 
f <- read.table(files[1], header = T)

# object that will hold the data
data <- data.frame(matrix(ncol = length(files), nrow = nrow(f)))
colnames(data) <- design_matrix$Description
rownames(data) <- f$Name
index = 1
# loop to grab data from NumReads column
for(file in files) {
  print(file)
  # read in data
  f <- read.table(file, header = T)
  # add to data
  data[,index] <- as.integer(f$NumReads)
  # advance index
  index = index +1
}
# go back to above directory

# make background gene list
rowsums <- rowSums(data)
background <- cbind(rownames(data), rowsums)
background <- background[which(background[,2] >= 10),]
fixed <- c()
for (id in background[,1]) {
  fix_ = unlist(strsplit(id, split = "\\."))[1]
  fixed <- c(fixed,fix_)
}
background <- cbind(background,fixed)
#write.csv(background, "background_genes.csv")


# make coldata factorial
rownames(design_matrix) = design_matrix$Description
design_matrix <- design_matrix[,3:6]
design_matrix$on.chip <- factor(design_matrix$on.chip)
design_matrix$flow <- factor(design_matrix$flow)
levels(design_matrix$flow) <- c("static", "flow")
design_matrix$huvec <- factor(design_matrix$huvec)
design_matrix$condition <- factor(design_matrix$condition)
# fix gene names and change to gene symbol
# quick function which takes the result matrix and transforms the names
gene_names <- function(results) {
  old <- rownames(results)
  fixed <- c()
  for (id in old) {
    fix_ = unlist(strsplit(id, split = "\\."))[1]
    fixed <- c(fixed,fix_)
  }
  
  new_names <- c()
  new_names <- mapIds(org.Hs.eg.db, 
                                    keys = fixed, 
                                    keytype="ENSEMBL", 
                                    column = "SYMBOL")
  
  for (entry in 1:length(new_names)) {
    if (is.na(new_names[entry])) {
      new_names[entry] <- "NA"
    }
  }
  
  
  rownames(results) <- new_names
  
  results <- results[-which(rownames(results) == "NA"),]
  
  output <- results
  return(output)
}

## Data selection ----
# to exclude the samples in wells
data <- data[,4:ncol(data)]
design_matrix <- design_matrix[4:nrow(design_matrix),]

# to exclude chip static samples 5,4 because they are on a different chip 
data <- data[,-c(5,4)]
design_matrix <- design_matrix[-c(5,4),]

background<- gene_names(background)

# Analysis -------------------------------------------------------------------
## DESeq2 differential expression setup ----
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = design_matrix,
                              design = ~ condition)

# filter out rowsums with less than 10 counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# this runs differential expression
dds <- DESeq(dds)

# Data transformations (different type of normalised data)
ntd <- normTransform(dds)  # just nromalisation
vsd <- vst(dds, blind = FALSE)  # variance stabilising transform
rld <- rlog(dds, blind = FALSE)  # log2 regularisation transform
rld <- gene_names(rld)

## 1 vs 4 - static+nohuvecs vs flow+huvecs ----
one_four <- results(dds, contrast=c("condition", "4", "1"))
one_four <- one_four[which(one_four$padj <= 0.05),]
one_four <- one_four[which(abs(one_four$log2FoldChange) > 1),]
one_four <- one_four[which(one_four$baseMean > 18),]
one_four <- gene_names(one_four)
#write.csv(one_four, "one_four.csv")

# Enrichment of GO terms
GO_res <- enrichGO(gene = rownames(one_four), 
                   OrgDb = "org.Hs.eg.db" , 
                   keyType = "SYMBOL" , ont = "BP", 
                   pAdjustMethod = "fdr", 
                   universe = rownames(background),
                   pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                   minGSSize = 30, maxGSSize = 500
                   )
fit <- plot(barplot(GO_res, showCategory = 15, font.size = 14, order=F))

test <- GO_res[4,]$geneID

test <- c(unlist(strsplit(test, split = "/")))

dds_counts <- counts(dds,normalized = T)
names <- colnames(dds_counts)
dds_counts <- gene_names(dds_counts)
dds_counts <- t(apply(dds_counts, 1, scale))
colnames(dds_counts) <- names
dds_counts <- cbind(dds_counts[,1:3], dds_counts[,10:12])

Heatmap(dds_counts[test,], cluster_rows = T, cluster_columns = T, 
        column_labels = colnames(dds_counts), name = "Z-score", 
        row_labels = test)

## 2 vs 4  - flow+nohuvecs vs flow+huvecs ----
two_four <- results(dds, contrast=c("condition", "4", "2"))
two_four <- two_four[which(two_four$padj <= 0.05),]
two_four <- two_four[which(abs(two_four$log2FoldChange) > 1),]
two_four <- two_four[which(two_four$baseMean > 18),]
two_four <- gene_names(two_four)
#write.csv(two_four, "two_four.csv")

# Enrichment of GO terms
GO_res <- enrichGO(gene = rownames(two_four), 
                   OrgDb = "org.Hs.eg.db" , 
                   keyType = "SYMBOL" , ont = "BP", 
                   pAdjustMethod = "fdr", 
                   universe = rownames(background),
                   pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                   #minGSSize = 30, maxGSSize = 500
)
fit <- plot(barplot(GO_res, showCategory = 15, font.size = 14, order=F))

## 3 vs 4 - static+huvecs vs flow+huvecs ----
three_four <- results(dds, contrast=c("condition", "4", "3"))
three_four <- three_four[which(three_four$padj <= 0.05),]
three_four <- three_four[which(abs(three_four$log2FoldChange) > 1),]
three_four <- three_four[which(three_four$baseMean > 18),]
three_four <- gene_names(three_four)
#write.csv(three_four, "three_four.csv")

# Enrichment of GO terms
GO_res <- enrichGO(gene = rownames(three_four), 
                   OrgDb = "org.Hs.eg.db" , 
                   keyType = "SYMBOL" , ont = "BP", 
                   pAdjustMethod = "fdr", 
                   universe = rownames(background),
                   pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                   minGSSize = 30, maxGSSize = 500
)
fit <- plot(barplot(GO_res, showCategory = 15, font.size = 14, order=F))

## 1 vs 2 - static+huvecs vs flow+huvecs ----
one_two <- results(dds, contrast=c("condition", "2", "1"))
one_two <- one_two[which(one_two$padj <= 0.05),]
one_two <- one_two[which(abs(one_two$log2FoldChange) > 1.5),]
one_two <- one_two[which(one_two$baseMean > 18),]
one_two <- gene_names(one_two)
#write.csv(one_two, "one_two.csv")

# Enrichment of GO terms
GO_res <- enrichGO(gene = rownames(one_two), 
                   OrgDb = "org.Hs.eg.db" , 
                   keyType = "SYMBOL" , ont = "BP", 
                   pAdjustMethod = "fdr", 
                   universe = rownames(background),
                   pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                   minGSSize = 30, maxGSSize = 500
)
fit <- plot(barplot(GO_res, showCategory = 15, font.size = 14, order=F))

## 1 vs 3 - static+huvecs vs flow+huvecs ----
one_three <- results(dds, contrast=c("condition", "3", "1"))
one_three <- one_three[which(one_three$padj <= 0.05),]
one_three <- one_three[which(abs(one_three$log2FoldChange) > 1.5),]
one_three <- one_three[which(one_three$baseMean > 18),]
one_three <- gene_names(one_three)
#write.csv(one_three, "one_three.csv")


GO_res <- enrichGO(gene = rownames(one_three), 
                   OrgDb = "org.Hs.eg.db" , 
                   keyType = "SYMBOL" , ont = "BP", 
                   pAdjustMethod = "fdr", 
                   universe = rownames(background),
                   pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                   minGSSize = 30, maxGSSize = 500)
fit <- plot(barplot(GO_res, showCategory = 15, font.size = 14, order=TRUE))

# side note for gideon : ZBTB16 is a gene for kidney development 
# and renal system development


## Other small comparisons ----
# we want a heatmap of genes from interesting GO terms with expression 
# from the lists 

cq_list <- read.csv("genes_list_1vs4.csv")
cq_list$genes <- mapIds(org.Hs.eg.db, 
                    keys = cq_list[,2], 
                    keytype="ENSEMBL", 
                    column = "SYMBOL")

# overlap between our two list
length(which(cq_list$genes %in% rownames(one_two)))

# positive regulation of blood vessel development
pos_blood <- c("PKM","HMOX1","CAV1","SERPINE1","COL1A1","SPARC","WNT5A","EGR1" ,
               "MMP19",  "EREG",  "ARHGAP22",  "CDH13",  "MMP14",  "ITGA5",  
               "ESM1",  "ANPEP",  "GREM1",  "CXCL8" , "HAS2","SPHK1")
length(which(pos_blood %in% rownames(one_two)))
length(which(pos_blood %in% cq_list$genes))

# Visualisation -----------------------------------------------------------
## Total dataset ----
# Heatmap of whole datatset using top 20 expressed genes

### PCA ----
vsdata <- varianceStabilizingTransformation(dds, blind=FALSE)
# rename the points to the correct names
vsdata@colData@listData[["condition"]] <- factor(c(rep("Static w/o vasc.", 3), 
                                                   rep("Flow w/o vasc.", 3), 
                                                   rep("Static with vasc.", 3), 
                                                   rep("Flow with vasc.", 3)))
plotPCA(vsdata, intgroup = "condition")


dists <- dist(t(assay(vsdata)))
plot(hclust(dists))

### Heatmap of all samples ----
# this is the LRT test, this is not what we use for the GO term analysis
# the nice thing about this is it gives you a p-value for genes across
# all conditions.
dds <- DESeq(dds,test = "LRT", reduced = ~1)
# this is top significant based on test LRT
sigs <- rowData(dds)
sigs <- sigs[which(sigs$LRTPvalue <= 0.01),]
top_100 <- assay(rld)[order(sigs$LRTPvalue, decreasing = F)[1:100],]

# try and get the top 20 fold change between all conditions which are sigs
cond_4_1 <- rownames(sigs[order(abs(sigs$condition_4_vs_1), decreasing = T),])[1:50]
cond_3_1 <- rownames(sigs[order(sigs$condition_3_vs_1, decreasing = T),])[1:25]
cond_2_1 <- rownames(sigs[order(sigs$condition_2_vs_1, decreasing = T),])[1:25]
all_conds <- c(cond_4_1,cond_3_1,cond_2_1)
top_100 <- assay(rld)[all_conds,]
top_100 <- t(scale(t(top_100)))
splits <- factor(c(rep(1:4, each = 3)), levels = c(1,3,2,4))


ha <- HeatmapAnnotation(
  foo1 = anno_block(gp = gpar(fill = 1:4),
                    labels = c("Static w/o vasc.", "Static with vasc.","Flow w/o vasc.","Flow with vasc." ), 
                    labels_gp = gpar(col = "white", fontsize = 8)))

# Column ordering
splits <- factor(c(rep(1:4, each = 3)), levels = c(1,3,2,4))
column_dend <- as.dendrogram(hclust(dist(t(top_100))))
column_dend <- rotate(column_dend, order = c(1:12))

# Row ordering
row_clusters <- kmeans(top_100,4)
row_splitting <- factor(paste0("Cluster\n", row_clusters$cluster), levels=c("Cluster\n4","Cluster\n3","Cluster\n2","Cluster\n1"))

Heatmap(top_100, name = "Z score",
        cluster_rows = T, 
        cluster_columns = T,
        #cluster_columns = column_dend, 
        #column_labels = colnames(top_100), 
        column_split = splits, 
        cluster_column_slices = F,
        cluster_row_slices = F,
        #column_order = colnames(top_100),
        #row_order = rownames(top_100),
        #row_labels = rownames(top_100),
        row_names_gp = gpar(fontsize = 3),
        row_split = row_splitting,
        #row_order = c(2,1,4,3),
        top_annotation = ha,
        border = T
        )

## One vs Four - Figures ---------------------------------------------------
names_one_four <- c("Static w/o vasc. 1","Static w/o vasc. 2", "Static w/o vasc. 3", 
           "Flow with vasc. 1","Flow with vasc. 2", "Flow with vasc. 3")
### Heatmap of top 20 significantly differentially expressed genes from rld data ----
#### lfcshrink ----
# take rlog data from condition 4 vs 1 and shrink it using lfcshrink function
resNorm <- lfcShrink(dds, type = "apeglm", coef = 4)
resNorm <- gene_names(resNorm)

# extract top 20 diff genes based on shrinkage
one_four_top_20 <- rownames(resNorm[order(abs(resNorm$log2FoldChange),decreasing = T),])[1:20]
# take count data from rlog transformed data
sig_resNorm <- assay(rld)[which(rownames(assay(rld)) %in% one_four_top_20),]
# subset to only have one and four
sig_resNorm <- sig_resNorm[,c(1:3,10:12)]
# scale data so that mean expression is set to 0
sig_resNorm <- t(scale(t(sig_resNorm)))

# Top annotation
ha <- HeatmapAnnotation(
  foo1 = anno_block(gp = gpar(fill = 1:4),
                    labels = c("Static w/o vasc.", "Flow with vasc." ), 
                    labels_gp = gpar(col = "white", fontsize = 8)))

# plot on heatmap
Heatmap(sig_resNorm, cluster_rows = T, cluster_columns = T, 
        column_labels = colnames(sig_resNorm), name = "Z score", 
        row_labels = rownames(sig_resNorm))

#### normal ----
# extract top 20 diff genes 
one_four_top_20 <- rownames(one_four[order(abs(one_four$log2FoldChange),decreasing = T),])[1:50]
# take count data from rlog transformed data
sig_resNorm <- assay(rld)[which(rownames(assay(rld)) %in% one_four_top_20),]
# subset to only have one and four
sig_resNorm <- sig_resNorm[,c(1:3,10:12)]
# scale data so that mean expression is set to 0
sig_resNorm <- t(scale(t(sig_resNorm)))

# Top annotation
ha <- HeatmapAnnotation(
  foo1 = anno_block(gp = gpar(fill = 1:2),
                    labels = c("Static w/o vasc.", "Flow with vasc." ), 
                    labels_gp = gpar(col = "white", fontsize = 8)))

column_splitting <- factor(c(rep("1", 3), rep("2", 3)), levels = c("1", "2"))
# plot heatmap
Heatmap(sig_resNorm, cluster_rows = T, cluster_columns = T,
        cluster_column_slices = F,
        column_labels = names_one_four,
        name = "Z score", 
        row_labels = rownames(sig_resNorm),
        top_annotation = ha,
        column_split = column_splitting
        )


### Heatmap of genes from GO terms (take counts from rld)
## Genes from blood circulation GO terms, extracted post factum and made into vector
blood_circ_term <- "ATP1A2,ATP1B1,ADORA1,AGT,BMP10,DES,KCNE4,HTR2B,P2RY1,SHOX2,CD38,CORIN,NPY1R,ADAMTS16,NPR3,EDN1,SGK1,CAV1,CACNB2,CXCL12,KCNMA1,SERPING1,OLR1,EDNRB,GCH1,C2CD4B,ANPEP,CACNA1H,TRPV1,ATP1B2,ACE,DSC2,DSG2,AZU1,FFAR3,CEACAM1,APOE,CTSZ,HMOX1,ADM2"
blood_circ_term <- unlist(strsplit(blood_circ_term, split = ","))

blood_term_resNorm <- assay(rld)[which(rownames(assay(rld)) %in% blood_circ_term),]
# subset to only have one and four
blood_term_resNorm <- blood_term_resNorm[,c(1:3,10:12)]
# subset to only the ones which are significant
which(blood_circ_term %in% rownames(one_four))
# scale data so that mean expression is set to 0
blood_term_resNorm <- t(scale(t(blood_term_resNorm)))
# change the colnames to fit the paper
colnames(blood_term_resNorm) <- names_one_four
# plot on heatmap
ha <- HeatmapAnnotation(
  foo1 = anno_block(gp = gpar(fill = 1:2),
  labels = c("Static", "Flow"), 
  labels_gp = gpar(col = "white", fontsize = 10)))

column_dend <- as.dendrogram(hclust(dist(t(blood_term_resNorm))))
column_dend <- rotate(column_dend, order = c(6:1))

Heatmap(blood_term_resNorm, cluster_rows = T, cluster_columns = column_dend, 
        column_labels = colnames(blood_term_resNorm),
        name = "Z score",
        row_labels = rownames(blood_term_resNorm),
        column_names_rot = 45,
        column_split = 2,top_annotation = ha)




### Exrtra cellular matrix organisation Go-term heatmap ----
extracell <- c("TNFRSF1B/PDPN/MATN1/COL9A2/CTSK/DPT/AGT/DPP4/MELTF/TLL1/ADAMTS16/ADAMTS12/FGFR4/ADAMTS2/TNXB/COL12A1/CAV1/TNFRSF11B/COL14A1/ADAMTSL1/COL15A1/NTNG2/ADAMTSL2/ITGA8/ADAMTS14/WT1/MMP1/COL2A1/MMP19/MMP17/MMP14/FBLN5/SPINT1/MFAP4/COL5A3/ADAMTS1/MMP11/FBLN1")
extracell <- unlist(strsplit(extracell, split = "/"))

extracell_rld <- assay(rld)[which(rownames(assay(rld)) %in% extracell),]
# subset to only have one and four
extracell_rld <- extracell_rld[,c(1:3,10:12)]
# scale data so that mean expression is set to 0
extracell_rld <- t(scale(t(extracell_rld)))
# change names to match paper
colnames(extracell_rld) <- names_one_four
# plot on heatmap
ha <- HeatmapAnnotation(
  foo1 = anno_block(gp = gpar(fill = 1:2),
                    labels = c("Static", "Flow"), 
                    labels_gp = gpar(col = "white", fontsize = 10)))
# rearrange cluster columns
column_dend <- as.dendrogram(hclust(dist(t(extracell_rld))))
column_dend <- rotate(column_dend, order = c(6:1))

Heatmap(extracell_rld, cluster_rows = T, cluster_columns = column_dend, 
        column_labels = colnames(extracell_rld), name = "Z score", 
        row_labels = rownames(extracell_rld),column_names_rot = 45,
        column_split = 2,top_annotation = ha)

### Extracellualr matrix organisation volcano plot ----
EnhancedVolcano(one_four, x = "log2FoldChange", y = "padj", lab = rownames(one_four),
                pCutoff = 0.05, 
                FCcutoff = 1,
                selectLab = extracell,
                max.overlaps = Inf,
                maxoverlapsConnectors = NULL,
                drawConnectors = TRUE,
                legendPosition = 'right')


### Heatmap of selected genes from GO terms ----
selected <- c("TIMP2",  "P4HA2",  "COL5A3",  "PLOD1",  "ADAMTS2",  "ITGA6",  
              "TIMP1",  "CAV1",  "SERPINE1",  "COL1A1",  "VWF",  "COL12A1",  
              "SPARC",  "COL7A1",  "MMP19",  "FBLN5",  "COL6A2",  "CTSK",  
              "ADAMTS12",  "MMP14",  "ITGA5",  "COL6A3",  "ITGA2",  "GREM1",  
              "HAS2",  "CD151",  "MMP1", "ADAMTS16", "KCNMA1", "CAV1", "ANPEP", 
              "HMOX1", "BMP10", "OLR1", "EDNRB", "AGT")

blood_circ_list <- c("KCNMA1 ADAMTS16 CAV1 ANPEP ACE HMOX1 SHOX2 BMP10 OLR1 EDNRB AGT")
blood_circ_list <- unlist(strsplit(blood_circ_list, split = " "))
extracell_list <- c("CTSK TLL1 ADAMTS12 ADAMTS2 COL12A1 ADAMTSL1 COL15A1 NTNG2 ADAMTS14 MMP1 MMP19 MMP17 MMP14 FBLN5 COL5A3 MMP11")
extracell_list <- unlist(strsplit(extracell_list, split = " "))
selected <- c(blood_circ_list,extracell_list)

# take out only these genes from one_four
best_genes <- assay(rld)[which(rownames(assay(rld)) %in% selected),]
# subset to only have one and four
best_genes <- best_genes[,c(1:3,10:12)]
# scale data so that mean expression is set to 0
best_genes <- t(scale(t(best_genes)))
# change names to match paper
colnames(best_genes) <- names_one_four

# This makes the top annotation which puts 2 blocks to label the columns 
# where samples are from
ha <- HeatmapAnnotation(
  foo1 = anno_block(gp = gpar(fill = 1:2),
  labels = c("Static", "Flow"), 
  labels_gp = gpar(col = "white", fontsize = 10)))
# rearrange cluster columns (this uses a package called dendextend)
column_dend <- as.dendrogram(hclust(dist(t(best_genes))))
column_dend <- rotate(column_dend, order = c(6:1))
splits <- factor(c(rep("Blood Circulation", length(blood_circ_list)), rep("Extracellular Matrix", length(extracell_list))), levels = c("Blood Circulation", "Extracellular Matrix"))
# this is the actual heatmap command
Heatmap(best_genes, cluster_rows = FALSE, cluster_columns = column_dend, 
        column_labels = colnames(best_genes), name = "Z score", 
        row_labels = rownames(best_genes),
        column_names_rot = 70,
        column_split = 2,top_annotation = ha, 
        row_dend_reorder = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = F,
        row_order = selected,
        #row_split = splits,
        border = F
        )


## this is based on the gene list pos_blood
blood_term_resNorm <- assay(rld)[which(rownames(assay(rld)) %in% pos_blood),]
# subset to only have one and four
blood_term_resNorm <- blood_term_resNorm[,c(1:3,10:12)]
# scale data so that mean expression is set to 0
blood_term_resNorm <- t(scale(t(blood_term_resNorm)))
# plot on heatmap
Heatmap(blood_term_resNorm, cluster_rows = T, cluster_columns = T, 
        column_labels = colnames(blood_term_resNorm), name = "Z score", 
        row_labels = rownames(blood_term_resNorm))




# Go term regulation of blood circulation
reg_blood <- GO_res[12,]$geneID



### Volcano Plot ----
EnhancedVolcano(one_four, x = "log2FoldChange", y = "padj", 
                lab = rownames(one_four),
                pCutoff = 0.05, 
                FCcutoff = 1,
                #selectLab = selected,
                #max.overlaps = Inf,
                #maxoverlapsConnectors = NULL,
                drawConnectors = TRUE,
                legendPosition = 'right',
                xlim = c(-10,10),
                title = "",
                subtitle = ""
                )



