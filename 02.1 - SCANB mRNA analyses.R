# SCAN-B cohort
# GSE81538

# Gene expression data in FPKM were generated using cufflinks 2.1.1 and pre-processed by collapsing on 
# 27,979 unique gene symbols (sum of FPKM values of each matching transcript), 
# keeping 18,802 NCBI RefSeq NM-category (mRNA) gene symbols and adding to each expression measurement 0.1 FPKM, 
# performing a log2 transformation.


library(tidyverse)
library(openxlsx)
library(conflicted)
conflict_prefer_all("dplyr")

setwd("C:/Users/idisb/Desktop/Andrés/02 - Projects/01 - TNBC HER2-low/SCAN-B")

data <- read_csv("GSE81538_gene_expression_405_transformed.csv")
metadata <- read.xlsx("Final Metadata TNBC HER2low GSE81538.xlsx")

myData <- data %>%
  column_to_rownames(var = "...1")
myData <- myData[,metadata$Sample.Name]
myData.norm <- 2^myData
myData.norm <- myData.norm - 0.1
myData.norm <- log2(myData.norm + 1) ### Ya tenemos los datos normalizados


metadata.filt <- metadata %>%
  mutate(Color = case_when(HER2_group %in% "HER2-zero" ~ "#884B99",
                           HER2_group %in% "HER2-low" ~ "#7CEBB0" )) %>%
  mutate(her2_clin.fin = her2_clin1 + her2_clin2 + her2_clin3) %>%
  mutate(HER2_group = case_when(her2_clin.fin == 0 ~ "HER2-zero",
                                her2_clin.fin >= 1 ~ "HER2-low"))

table(metadata.filt$HER2_group)

RNA.anot <- read_tsv("Gene_metadata_annotations.txt")  %>%
  filter(gene_name %in% rownames(myData.norm)) %>%
  select(gene_name, gene_id, gene_type)

myData.norm <- myData.norm[, metadata.filt$Sample.Name]
myData.HER2low <- myData.norm[,metadata.filt$Sample.Name[metadata.filt$HER2_group == "HER2-low"]] 
myData.HER2zero <- myData.norm[,metadata.filt$Sample.Name[metadata.filt$HER2_group == "HER2-zero"]]


# Normality test
shapiro.test(as.numeric(myData.norm[1,]))

# We will use a t-test to calculate statistical differences, then adjust for FDR

# Log2FC will be calculated as mean(log2(condition1)) - mean(log2(condition2))
# The data is already in log2 format https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5826329

# Calculate the mean of each gene
HER2low.mean <- apply(myData.HER2low, 1, mean)
HER2zero.mean <- apply(myData.HER2zero, 1, mean)
baseMean <- apply(myData.norm, 1, mean)

log2FC <- HER2low.mean - HER2zero.mean

# Compute statistical significance #
pvalue = NULL
tstat = NULL
for(i in 1 : nrow(myData.HER2low)) {
  x = myData.HER2low[i,]
  y = myData.HER2zero[i,]
  
  result <- tryCatch({
    t = t.test(as.numeric(x), as.numeric(y))
    list(p.value = t$p.value, statistic = t$statistic)
  }, error = function(e) {
    list(p.value = NA, statistic = NA) # Devuelve NA si hay un error
  })
  
  pvalue[i] = result$p.value
  tstat[i] = result$statistic
}

p.adj <- p.adjust(pvalue, method = "fdr", n = length(pvalue)) # Correction by False Discovery Rate

res.naive <- as.data.frame(cbind(baseMean, pvalue, p.adj, log2FC)) %>%
  rownames_to_column(var = "gene_name") %>%
  drop_na() %>%
  left_join(RNA.anot, by = "gene_name") %>%
  filter(!duplicated(gene_name))



# GSEA
# Enrichment analyses #
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

# 2. WITH ALL THE GENES --> We need the background: THIS ONE!!
# To make the list: we want the log2 fold change 
GSEA.list.df <- res.naive %>%
  drop_na(gene_id)

GSEA.list <- GSEA.list.df$log2FC
names(GSEA.list) <- GSEA.list.df$gene_id

# omit any NA values 
GSEA.list <- na.omit(GSEA.list)

# sort the list in decreasing order (required for clusterProfiler)
GSEA.list <- sort(GSEA.list, decreasing = TRUE)

#We eliminate ... (variants of genes)
names(GSEA.list) <- sub("\\..*", "", names(GSEA.list))

# Gene Set Enrichment Analysis of Gene Ontology
set.seed(231)
GSEA.GO <- gseGO(geneList=GSEA.list, 
                 ont ="BP", 
                 keyType = "ENSEMBL", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = "org.Hs.eg.db",
                 nPermSimple = 10000,
                 pAdjustMethod = "BH",
                 seed = TRUE,
                 eps = 0)

GSEA.GO.df <- as.data.frame(GSEA.GO)

## Plots
gse.paper <- GSEA.GO.df[1:10,]


gse.paper <- gse.paper %>%
  mutate(Description = fct_reorder(Description, desc(NES))) %>%
  arrange(NES)

ggplot(gse.paper, aes(y = Description, x = NES)) +
  geom_segment(aes(yend = Description, xend = 0)) +
  geom_point(size = 6, color ="black", fill = alpha("orange", 0.5), shape = 21) + 
  theme_bw() + 
  labs(x = "NES (Normalized Enrichment Score)")

which(GSEA.GO$Description == "lymphocyte activation")

gseaplot(GSEA.GO,geneSetID=102,title=GSEA.GO$Description[102]) #GSEA



## IMMUNE-RELATED
gse.paper <- GSEA.GO.df %>%
  filter(Description %in% c("immune system development",
                            "B cell activation",
                            "mitotic cell cycle checkpoint signaling",
                            "B cell differentiation",
                            "antigen receptor-mediated signaling pathway",
                            "lymphocyte activation",
                            "regulation of leukocyte cell-cell adhesion",
                            "positive regulation of lymphocyte activation",
                            "regulation of lymphocyte activation"))

gse.paper <- gse.paper %>%
  mutate(Description = fct_reorder(Description, desc(NES))) %>%
  arrange(NES)

ggplot(gse.paper, aes(y = Description, x = NES)) +
  geom_segment(aes(yend = Description, xend = 0)) +
  geom_point(size = 6, color ="black", fill = alpha("orange", 0.5), shape = 21) + 
  theme_bw() + 
  labs(x = "NES (Normalized Enrichment Score)")



# HLAs cascade plot
HLAs2 <- res.naive %>%
  filter(grepl("HLA-", gene_name)) %>%
  filter(gene_type %in% "protein_coding") %>%
  mutate(sig = case_when(pvalue < 0.05  ~ "Y",
                         T ~ "N"))

# Crec que queda prou bé, es colors havia pensat posarlos a nes peu de pàgina
# directament

ggplot(HLAs2, aes(x= fct_reorder(gene_name,desc(log2FC)), y = log2FC, fill = sig)) +
  geom_col(color = "black", width = 1) + 
  ylim(-0.6, 0.6) +
  theme_classic() + 
  scale_fill_manual(values = c("N" = "gray", "Y" = "tomato")) + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1))




