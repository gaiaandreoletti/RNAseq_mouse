##Gaia Andreoletti
## Dec 2021

setwd("/zebrafish_data/")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("tximport")
library("tximport")
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("/s3/ica/folder/Danio_rerio.GRCz11.109.gtf")
keytypes(txdb)
columns(txdb)
k <- keys(txdb, keytype = "GENEID")
#AnnotationDbi::select(txdb, keys = k, columns="GENEID", keytype="TXNAME")
df <- AnnotationDbi::select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1]
write.csv(tx2gene, "tx2gene_GRCz11.109.csv")
head(tx2gene)
tx2gene[grep("149835", tx2gene$TXNAME), ]

# ## mtm1 zebrafish trascripts 
# TXNAME             GENEID
# ENSDART00000149835 ENSDARG00000037560
# BiocManager::install("org.Dr.eg.db") #genome wide annotaion for Zebrafish more comprehensive
library(org.Dr.eg.db)
keytypes(org.Dr.eg.db)
# ens <-  rownames
mapIds(org.Dr.eg.db, keys = k,
       column = c('SYMBOL'), keytype = 'ENSEMBL')
k <- keys(org.Dr.eg.db, keytype = "ENSEMBLTRANS")
df <- AnnotationDbi::select(org.Dr.eg.db, keys = k,
       columns = c('ENSEMBLTRANS',  'ENTREZID', 'SYMBOL', 'ENSEMBLTRANS'),
       keytype = 'ENSEMBLTRANS')
head(df)
tx2gene <- df[,cbind(1,3)]
head(tx2gene)
dim(tx2gene)
tx2gene <- na.omit(tx2gene)
grep("mtm1", tx2gene$SYMBOL)
tx2gene[grep("^mtm1", tx2gene$SYMBOL), ]

# BiocManager::install("AnnotationHub")
# library(AnnotationHub)
# hub <- AnnotationHub()
# query(hub, c("Danio rerio","txdb"))
# ssdb <- hub[["AH52254"]]
# ssdb
# genes(ssdb)

# a<- gsub("\\..*","",tx2gene[,1])
# b<- gsub("\\..*","",tx2gene[,2])
# c<-cbind(a,b)
# colnames(c)=colnames(tx2gene)
# tx2gene <- as.data.frame(c)
# head(tx2gene)

### abcb11b_d6 zebrafish ######################
dir= getwd()
samples <- read.table("Sample_info.txt", header = T)
files <- file.path(dir, "quants", samples$fastq_ID., "quant.sf")
# names(files) <- paste0("sample", 1:6)
names(files) <- samples$Sample_ID
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
# txi.salmon <- tximport(files, type = "salmon")
head(txi.salmon$counts)

####  asked for a PA of old cohort and BSEP data 6 march 2023 ############
write.csv(txi.salmon$counts, file ="BSEP_counts.csv", quote = FALSE)
old_cohort <- read.csv("26_Apr_2021_counts_zebrafish_no_zpr1.csv", header = TRUE, row.names = 1)
BSEP_counts <- txi.salmon$counts
tails_counts <- read.csv("../zebrafish_quantsec/Zebrafish_quantseq_countsMatrix.csv", header = TRUE, row.names = 1)
merged_counts <- merge(as.data.frame(old_cohort), as.data.frame(BSEP_counts), by="row.names", sort=FALSE)
head(merged_counts)
rownames(merged_counts) <- merged_counts$Row.names
merged_counts$Row.names <- NULL
head(merged_counts)
samples <- read.table("Sample_info_merged.txt", header = T)

library(DESeq2)
data_counts <- round(merged_counts,2)
ncol(data_counts) == nrow(samples)
rownames(samples) <-samples$Sample_ID
data_counts <- data_counts[,order(colnames(data_counts))]
samples <- samples[order(rownames(samples)),]
colnames(data_counts) == rownames(samples)

samples$condition <- as.factor(samples$condition)

dds <- DESeqDataSetFromMatrix(countData = round(data_counts),
                              colData = samples,
                              design = ~  condition)


cts = counts(dds)
geoMeans = apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds = estimateSizeFactors(dds, geoMeans=geoMeans)
dds <- dds[ rowSums( counts(dds) ) > 0 , ]
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE, fitType='local')
DESeq2::plotPCA(vsd, intgroup="condition")+ ggtitle(" PCA of WT, EtOH, Overfed, MTM and BSEP") 
ggsave(" PCA_WT_EtOH_Overfed_MTM_BSEP.pdf", width = 20, height = 20, units = "cm")

######################END #####################################
# BiocManager::install("DESeq2")
library(DESeq2)
data_counts <- round(txi.salmon$counts,2)
ncol(data_counts) == nrow(samples)
rownames(samples) <-samples$Sample_ID
data_counts <- data_counts[,order(colnames(data_counts))]
samples <- samples[order(rownames(samples)),]
colnames(data_counts) == rownames(samples)

samples$WT_KO <- as.factor(samples$WT_KO)

dds <- DESeqDataSetFromMatrix(countData = round(data_counts),
                              colData = samples,
                              design = ~  WT_KO)


cts = counts(dds)
geoMeans = apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds = estimateSizeFactors(dds, geoMeans=geoMeans)
dds <- dds[ rowSums( counts(dds) ) > 0 , ]
dds <- DESeq(dds)
summary(dds)

## plot mtm1 exp all samples #ENSDART00000149835
d <- plotCounts(dds, gene="mtm1", intgroup="WT_KO", returnData=TRUE, normalized = TRUE)
ggplot(d, aes(WT_KO, count, color = WT_KO)) + geom_point()
d <- plotCounts(dds, gene="abcb11b", intgroup="WT_KO", returnData=TRUE, normalized = TRUE)
ggplot(d, aes(WT_KO, count, color = WT_KO)) + geom_point()

d$Condition <- factor(d$Condition, levels = c("MtDMSO1", "MtDMSO2", "MtDMSO3" ,"WTDMSO1", "WTDMSO2", "WTDMSO3"))
ggboxplot(d, x = "WT_KO", y = "count", 
          color = "WT_KO", palette =c("#00AFBB", "#E7B800"),
          add = "jitter", shape = "WT_KO") + ggtitle("abcb11b expression across all samples") + theme_bw() +
  stat_compare_means(label.y = 499)                   # Add global p-value
ggsave("mouse_cohort_abcb11b_exp.pdf", height=6, width=5.5)


# dds <- estimateSizeFactors(dds)
# nsamp <- 17 #80% of 22 
# keep <- rowSums(counts(dds, normalized=TRUE)> 5)>= nsamp 
# table(keep) #  11800   keep (7332  discarded) (( for 21 samples 7437 false, 11681 keep))
# dds <- dds[keep,]  
# dds <- DESeq(dds, fitType='local')
# summary(dds)
# 
####  "mtm1_mutants","WT" ####
resultsNames(dds) 
res1 <- results(dds, contrast = c("WT_KO","WT","Mut"))
head(res1,2)
dim(res1)
table(res1$padj<0.05)#1234      
res1 <- res1[order(res1$padj), ] 
head(res1,2)
res1_data <- merge(as.data.frame(res1), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(res1_data)[1] <- "Gene" 
head(res1_data,2) 
table(res1$padj<0.05 & res1$log2FoldChange > 1)# 433 in WT
table(res1$padj<0.05 & res1$log2FoldChange < - 1)# 599 in KO
library(biobroom)
library(dplyr)
write.csv(res1, file=paste(date_is,"diffexpr-results_ zebrafish_data_wt_vs_mut.csv"))
vsd <- vst(dds, blind=FALSE, fitType='local')
DESeq2::plotPCA(vsd, intgroup="WT_KO")+ ggtitle("WT_KO") 
 
# filter(df, padj < 0.05 & log2FoldChange > 1)
res1_data %>% 
  dplyr::filter(padj < 0.05 & log2FoldChange> 1) %>% 
  dplyr::select(Gene) %>% 
  knitr::kable(.)

res1_data %>% 
  dplyr::filter(padj < 0.05 & log2FoldChange < -1) %>% 
  dplyr::select(Gene) %>% 
  knitr::kable(.)

# DESeq2::plotPCA(vsd, intgroup="sample_name")+ ggtitle("sample_name")
# DESeq2::plotPCA(vsd, intgroup="Condition")+ ggtitle("Condition")
# DESeq2::plotPCA(vsd, intgroup="Batch")+ ggtitle("Batch")
# # assay(vsd)<- limma::removeBatchEffect(assay(vsd), vsd$Lane)
# dev.off()
library(pheatmap)
pdf("zebrafish_data_heatmap.pdf")
padj.cutoff <- 0.05 
lfc.cutoff <- 1 
threshold <- res1$padj < padj.cutoff & abs(res1$log2FoldChange)> lfc.cutoff 
length(which(threshold))
res1$threshold <- threshold     
sigOE <- data.frame(subset(res1, threshold==TRUE))
normalized_counts <- counts(dds, normalized=T)
norm_OEsig <- normalized_counts[sigOE$gene,] 
select <- rownames(subset(res1, threshold==TRUE))
df <- as.data.frame(colData(dds)[,"WT_KO"])
rownames(df)<- colData(dds)$Sample_ID
#pheatmap(assay(vsd)[select,], show_rownames=T, annotation_col=df,cutree_cols = 2, cutree_rows = 2, fontsize = 4, fontsize_row  = 6, scale = "none",width = 2, height = 1)
names <- c("MTM1_R1_resub_Jan11_21_S20","MTM1_set_3_S13",
           "MTM1_3R_resub_Jan18_21_S22",  "MTM1_set_1_S11","MTM1_set_2_S12", 
           "WT_EF_set_1_S5", "WT_EF_set_2_S6", "WT_EF_set_3_S7","WT_EF_set_4_S8")
library(pheatmap)
pheatmap(assay(vsd)[select,], show_rownames=F, annotation_col=df, fontsize = 4, fontsize_row  = 6, scale = "row",width = 2, height = 1)
#pheatmap(assay(vsd)[select,names], show_rownames=F, annotation_col=df, fontsize = 4, fontsize_row  = 6, scale = "row",width = 2, height = 1)
dev.off()

# d<-assay(vsd)[select,names]
# write.table(d,"table.txt", quote = FALSE)
# d <- read.table("table.txt", header = T, row.names = 1)
# pheatmap(d, show_rownames=T,  fontsize = 4, fontsize_row  = 6, scale = "none",width = 2, height = 1)

#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

pdf("zebrafish_data_heatmap.pdf_Volcano_padj0.05_log2FC2.pdf")	  
EnhancedVolcano(res1, lab = rownames(res1), 
                x = 'log2FoldChange',      
                y = 'padj', 
                xlim = c(-20, 20),     
                pCutoff = 0.05, 
                #selectLab = names_uniq$V1, 
                FCcutoff = 2)
                
dev.off()

# DESeq results to pathways in 60 Seconds with the fgsea package
library(fgsea)
ens2symbol = mapIds(org.Dr.eg.db, 
                    keys=rownames(res1),  
                    column="ENSEMBL", 
                    keytype="SYMBOL",  
                    multiVals="first")

ranks <- res1_data$log2FoldChange
names(ranks) <- res1_data$Gene
head(ranks, 20)

barplot(sort(ranks, decreasing = T))

## https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
# install.packages("msigdbr")
library(msigdbr)
msigdbr_species()
all_gene_sets = msigdbr(species = "Danio rerio")
head(all_gene_sets)
msigdbr_collections()

all_gene_sets %>%
  dplyr::filter(gs_cat == "H") %>%
  head()

# pathways.hallmark <- gmtPathways("h.all.v7.4.symbols.gmt")
# pathways.hallmark %>% 
#   head() %>% 
#   lapply(head)

h_gene_sets = msigdbr(species = "Danio rerio", category = "H")
head(h_gene_sets)

msigdbr_list = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)
fgseaRes <- fgsea(pathways=msigdbr_list, stats=ranks, nperm=1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="up/down reg pathways zebrafish_data WT vs Mut") + 
  theme_minimal()

topUp <- fgseaRes %>% 
  dplyr::filter(ES > 0) %>% 
  top_n(10, wt=-padj)

topDown <- fgseaRes %>% 
  dplyr::filter(ES < 0) %>% 
  top_n(10, wt=-padj)

topPathways <- bind_rows(topUp, topDown) %>% 
  arrange(-ES)

library(data.table)
fwrite(topPathways, file ="pathways zebrafish_data WT vs Mut.csv")

plotGseaTable(pathways.hallmark[topPathways$pathway], 
              ranks, 
              fgseaRes, 
              gseaParam = 0.5)


#### "Quantseq data" ####
setwd("/fast/analysis/20210406_project__Zebrafish/DataForRebuttal/zebrafish_quantsec")
dir= getwd()
samples <- read.table("Sample_info.txt", header = T)

dir= getwd()
samples <- read.table("Sample_info.txt", header = T)
count_matrix <- read.csv("Zebrafish_quantseq_countsMatrix.csv", header=T, row.names=1)
head(count_matrix)
dim(count_matrix)

# BiocManager::install("DESeq2")
library(DESeq2)
data_counts <- round(count_matrix,2)
ncol(data_counts) == nrow(samples)
rownames(samples) <-samples$Sample_ID
colnames(data_counts) <- c("MTM_KO_Tails1","MTM_KO_Tails2", "MTM_KO_Tails3", "MTM_KO_Tails4", "WT_Tails1","WT_Tails2","WT_Tails3" ,"WT_Tails4", "WT_Tails5")
data_counts <- na.omit(data_counts)
data_counts <- data_counts[,order(colnames(data_counts))]
samples <- samples[order(rownames(samples)),]
colnames(data_counts) == rownames(samples)
samples$WT_KO <- as.factor(samples$WT_KO)

dds <- DESeqDataSetFromMatrix(countData = round(data_counts),
                              colData = samples,
                              design = ~  WT_KO)


cts = counts(dds)
geoMeans = apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds = estimateSizeFactors(dds, geoMeans=geoMeans)
dds <- dds[ rowSums( counts(dds) ) > 0 , ]
dds <- DESeq(dds)
summary(dds)

## plot mtm1 exp all samples #ENSDART00000149835
# ## mtm1 zebrafish trascripts 
# TXNAME             GENEID
# ENSDART00000149835 ENSDARG00000037560

d <- plotCounts(dds, gene="ENSDART00000149835", intgroup="WT_KO", returnData=TRUE, normalized = TRUE)
ggplot(d, aes(WT_KO, count, color = WT_KO)) + geom_point()
# d <- plotCounts(dds, gene="abcb11b", intgroup="WT_KO", returnData=TRUE, normalized = TRUE)
# ggplot(d, aes(WT_KO, count, color = WT_KO)) + geom_point()

# d$Condition <- factor(d$Condition, levels = c("MtDMSO1", "MtDMSO2", "MtDMSO3" ,"WTDMSO1", "WTDMSO2", "WTDMSO3"))
ggboxplot(d, x = "WT_KO", y = "count", 
          color = "WT_KO", palette =c("#00AFBB", "#E7B800"),
          add = "jitter", shape = "WT_KO") + ggtitle("mtm1 expression across quantseq samples") + theme_bw() +
  stat_compare_means(label.y = 100)                   # Add global p-value
ggsave("quantseq_mtm1_exp.pdf", height=6, width=5.5)

#### "WT", "WT_KO" ####
resultsNames(dds) 
res1 <- results(dds, contrast = c("WT_KO","WT","MTM_KO"))
head(res1,2)
dim(res1)
table(res1$padj<0.05)# 8627     
res1 <- res1[order(res1$padj), ] 
head(res1,2)
res1_data <- merge(as.data.frame(res1), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(res1_data)[1] <- "ensGene" 
head(res1_data,2) 
res1_data_annotated <- res1_data %>% inner_join(tx2gene, by=c("ensGene"="TXNAME")) 
head(res1_data_annotated,2) 
table(res1$padj<0.05 & res1$log2FoldChange > 1)# 8558 in WT
table(res1$padj<0.05 & res1$log2FoldChange < - 1)# 69 in KO
library(biobroom)
library(dplyr)
write.csv(res1, file=paste(date_is,"diffexpr-results_zebrafish_quantsec_wt_vs_mut.csv"))
write.csv(res1_data_annotated, file=paste(date_is,"diffexpr-results_zebrafish_quantsec_wt_vs_mut_annotated.csv"))

vsd <- vst(dds, blind=FALSE, fitType='local')
DESeq2::plotPCA(vsd, intgroup="WT_KO")+ ggtitle("WT_KO") 
ggsave(" zebrafish_quantsec_rnaseq_fastq_zebrafish_cohort_pca.pdf", height=6, width=5.5)

grep("mtm1", res1_data_annotated$SYMBOL)
res1_data_annotated[grep("mtm1", res1_data_annotated$SYMBOL), ]

# filter(df, padj < 0.05 & log2FoldChange > 1)
res1_data_annotated %>% 
  dplyr::filter(padj < 0.05 & log2FoldChange> 1) %>% 
  dplyr::select(SYMBOL) %>% 
  knitr::kable(.)

res1_data_annotated %>% 
  dplyr::filter(padj < 0.05 & log2FoldChange < -1) %>% 
  dplyr::select(SYMBOL) %>% 
  knitr::kable(.)

# DESeq2::plotPCA(vsd, intgroup="sample_name")+ ggtitle("sample_name")
# DESeq2::plotPCA(vsd, intgroup="Condition")+ ggtitle("Condition")
# DESeq2::plotPCA(vsd, intgroup="Batch")+ ggtitle("Batch")
# # assay(vsd)<- limma::removeBatchEffect(assay(vsd), vsd$Lane)
# dev.off()
library(pheatmap)
pdf("quantseq_rnaseq_fastq_zebrafish_heatmap.pdf")
padj.cutoff <- 0.05
lfc.cutoff <- 1
threshold <- res1$padj < padj.cutoff & abs(res1$log2FoldChange)> lfc.cutoff 
length(which(threshold))
res1$threshold <- threshold     
sigOE <- data.frame(subset(res1, threshold==TRUE))
normalized_counts <- counts(dds, normalized=T)
norm_OEsig <- normalized_counts[sigOE$gene,] 
select <- rownames(subset(res1, threshold==TRUE))
df <- as.data.frame(colData(dds)[,"WT_KO"])
rownames(df)<- colData(dds)$Sample_ID
#pheatmap(assay(vsd)[select,], show_rownames=T, annotation_col=df,cutree_cols = 2, cutree_rows = 2, fontsize = 4, fontsize_row  = 6, scale = "none",width = 2, height = 1)
library(pheatmap)
pheatmap(assay(vsd)[select,], show_rownames=F, annotation_col=df, fontsize = 4, fontsize_row  = 6, scale = "row",width = 2, height = 1)
#pheatmap(assay(vsd)[select,names], show_rownames=F, annotation_col=df, fontsize = 4, fontsize_row  = 6, scale = "row",width = 2, height = 1)
dev.off()

# d<-assay(vsd)[select,names]
# write.table(d,"table.txt", quote = FALSE)
# d <- read.table("table.txt", header = T, row.names = 1)
# pheatmap(d, show_rownames=T,  fontsize = 4, fontsize_row  = 6, scale = "none",width = 2, height = 1)

#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

pdf("quantseq_zebrafish_heatmap.pdf_Volcano_padj0.05_log2FC2.pdf")	  
EnhancedVolcano(res1_data_annotated, 
                lab = res1_data_annotated$SYMBOL, 
                x = 'log2FoldChange',      
                y = 'padj', 
                xlim = c(-20, 20),     
                pCutoff = 0.05, 
                # selectLab = res1_data_annotated$SYMBOL, 
                FCcutoff = 2)

dev.off()

# DESeq results to pathways in 60 Seconds with the fgsea package
head(res1_data_annotated)
library(fgsea)
ens2symbol = mapIds(org.Dr.eg.db, 
                    keys=rownames(res1),  
                    column="ENSEMBL", 
                    keytype="ENSEMBLTRANS",  
                    multiVals="first")
ens2symbol <- as.data.frame(ens2symbol)

ranks <- res1_data_annotated$log2FoldChange
names(ranks) <- res1_data_annotated$GENEID
head(ranks, 20)

barplot(sort(ranks, decreasing = T))

## https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
# install.packages("msigdbr")
library(msigdbr)
msigdbr_species()
all_gene_sets = msigdbr(species = "Danio rerio")
head(all_gene_sets)
msigdbr_collections()

all_gene_sets %>%
  dplyr::filter(gs_cat == "H") %>%
  head()

# pathways.hallmark <- gmtPathways("h.all.v7.4.symbols.gmt")
# pathways.hallmark %>% 
#   head() %>% 
#   lapply(head)

h_gene_sets = msigdbr(species = "Danio rerio", category = "H")
head(h_gene_sets)

msigdbr_list = split(x = h_gene_sets$ensembl_gene, f = h_gene_sets$gs_name)
fgseaRes <- fgsea(pathways=msigdbr_list, stats=ranks, nperm=1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="up & down reg pathways quantseq fastq_zebrafish WT vs Mut") + 
  theme_minimal()

topUp <- fgseaRes %>% 
  dplyr::filter(ES > 0) %>% 
  top_n(10, wt=-padj)

topDown <- fgseaRes %>% 
  dplyr::filter(ES < 0) %>% 
  top_n(10, wt=-padj)

topPathways <- bind_rows(topUp, topDown) %>% 
  arrange(-ES)

library(data.table)
fwrite(topPathways, file ="pathways quantseq_rnaseq_fastq_zebrafish WT vs Mut.csv")

plotGseaTable(pathways.hallmark[topPathways$pathway], 
              ranks, 
              fgseaRes, 
              gseaParam = 0.5)

# ######found online https://rpubs.com/chapmandu2/171024#####
# setwd("/fast/analysis/20210406_project__Zebrafish/DataForRebuttal/zebrafish_quantsec")
# fastq_dir <- '/fast/analysis/20210406_project__Zebrafish/DataForRebuttal/zebrafish_quantsec/raw_fastq/'
# ref_fa <- 'Danio_rerio/Ensembl/Zv9/Sequence/WholeGenomeFasta/genome.fa'
# ref_dir <- 'Danio_rerio/Ensembl/Zv9/Sequence/WholeGenomeFasta/'
# gtf.file <- 'Danio_rerio/Ensembl/Danio_rerio.Zv9.77.gtf.gz'
# ???sqlite_file <- 'Homo_sapiens.GRCh38.84.sqlite'
# sqlite_path <- file.path(ref_dir, sqlite_file)
# workflow_stats <- list()
# 
# qa_plots <- function(shortfq_obj) {
#   x <- ShortRead::qa(shortfq_obj, lane='')
#   nuc_plot <- ggplot(x[['perCycle']]$baseCall, aes(x=Cycle, y=Count, colour=Base)) + geom_line() + cowplot::theme_cowplot()
#   qual_plot <- x[['perCycle']]$quality %>% group_by(Cycle) %>% summarise(QualScore=sum(Score*Count)/sum(Count)) %>% ungroup() %>% 
#     ggplot(aes(x=Cycle, y=QualScore)) + geom_line() + cowplot::theme_cowplot()
#   size_plot <- ggplot(data.frame(frag_size=width(shortfq_obj)), aes(x=frag_size)) + 
#     geom_bar(fill='pink', colour='black') + cowplot::theme_cowplot()
#   base_plot <- ggplot(data.frame(idx=1:50, dat=as.character(sread(shortfq_obj)[1:50]))) + 
#     geom_text(aes(x=1, y=idx, label=dat), size=rel(2.5), family='mono', hjust=0) + xlim(1,1.8) +
#     theme_bw() + theme(axis.text.x=element_text(size=0)) + xlab('')
#   return(list(nuc_plot=nuc_plot, qual_plot=qual_plot, size_plot=size_plot, base_plot=base_plot))
# } 
# 
# library("ShortRead")
# fq_fn <- grep('Zebrafish_', list.files(fastq_dir, 'fastq.gz'), value = TRUE)[1]
# fq_path <- file.path(fastq_dir, fq_fn)
# fq_data <- ShortRead::yield(FastqSampler(fq_path, 1000000))
# qp <- qa_plots(fq_data)
# cowplot::plot_grid(qp$nuc_plot, qp$qual_plot, qp$size_plot, labels='auto')
# qp$base_plot
# ShortRead:::.freqSequences(qa(fq_data, lane=''), "read")
# 
# fq <- fq_data
# fq <- ShortRead::narrow(fq, start=13)
# workflow_stats$step01_start <- length(fq)
# length(fq)
# 
# qp <- qa_plots(fq)
# cowplot::plot_grid(qp$nuc_plot, qp$qual_plot, qp$size_plot, labels='auto')
# qp$base_plot
# ShortRead:::.freqSequences(qa(fq_data, lane=''), "read")
# 
# fq <- ShortRead::trimTailw(fq, k=6, a="4", halfwidth=6)
# fq <- fq[width(fq) >= 32]
# workflow_stats$step02_poorqual <- length(fq)
# length(fq)
# 
# qp <- qa_plots(fq)
# cowplot::plot_grid(qp$nuc_plot, qp$qual_plot, qp$size_plot, labels='auto')
# 
# narrow_ranges <- ShortRead::trimEnds(sread(fq), right=TRUE, "A", relation="==", ranges=TRUE) 
# fq <- ShortRead::narrow(fq, start(narrow_ranges), end(narrow_ranges))
# fq <- fq[width(fq) >= 32]
# workflow_stats$step03_polya <- length(fq)
# length(fq)
# 
# qp <- qa_plots(fq)
# cowplot::plot_grid(qp$nuc_plot, qp$qual_plot, qp$size_plot, labels='auto')
# 
# polyApos <- regexpr(pattern= 'AAAAAAAAAA', sread(fq), fixed=TRUE) %>% unlist()
# polyAclip_idx <- which(polyApos >= 0)
# polyAclip <- width(fq)
# polyAclip[polyAclip_idx] <- polyApos[polyAclip_idx] 
# fq <- narrow(fq, end=polyAclip)
# fq <- fq[width(fq) >= 32]
# workflow_stats$step04_3padapt <- length(fq)
# length(fq)
# 
# qp <- qa_plots(fq)
# cowplot::plot_grid(qp$nuc_plot, qp$qual_plot, qp$size_plot, labels='auto')
# 
# library(RSubread)
# tmp_fastq <- tempfile()
# writeFastq(fq, tmp_fastq)
# tmp_bam <- tempfile()
# align(index=file.path(ref_dir,"zebrafish_subread_index"),
#       readfile1=tmp_fastq,   
#       input_format='gzFASTQ', 
#       output_file=tmp_bam,   
#       nthreads=8)
