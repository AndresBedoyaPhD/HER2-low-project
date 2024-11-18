library(openxlsx)
library(tidyverse)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(conflicted)
conflict_prefer_all("dplyr")

setwd("C:/Users/idisb/Desktop/Andrés/02 - Projects/01 - TNBC HER2-low")

# Upload Metadata
metadata <- read.xlsx("TNBC_TCGA_Final.xlsx")
metadata$Sample.ID <- paste0(metadata$Sample.ID, "A")
metadata <- metadata %>%
           mutate(Group = case_when(
             BC_Subtype_HER2 == "4.2_TNBC_HER2low" ~ "HER2-low",
             BC_Subtype_HER2 =="4.1_TNBC_HER2_0" ~ "HER2-zero"))


meta_HER2low <- metadata %>%
  filter(Group %in% "HER2-low")
meta_HER2zero <- metadata %>%
  filter(Group %in% "HER2-zero") 

# DNA METHYLATION - 450K #
# Search TCGA data filtered by barcode = Sample ID of curated metadata. 
query_DNA_methyl <- GDCquery(project = "TCGA-BRCA",
                             data.category = "DNA Methylation",
                             platform = "Illumina Human Methylation 450",
                             access = "open",
                             data.type = "Masked Intensities",
                             barcode = metadata$Sample.ID,
                             sample.type = "Primary Tumor")

data_DNAm <- getResults(query_DNA_methyl)

# Download data 
GDCdownload(query_DNA_methyl)

# Prepare and Normalize data by sesame
library(sesame)
dna.meth <- GDCprepare(query_DNA_methyl, summarizedExperiment = TRUE) 

# extract gene and sample metadata from summarizedExperiment object
coldata <- as.data.frame (colData(dna.meth))

# add condition: HER2low vs Her2-zero
coldata <- coldata %>%
  mutate(condition = case_when(coldata$sample %in% meta_HER2low$Sample.ID ~ "HER2-low",
                               coldata$sample %in% meta_HER2zero$Sample.ID ~ "HER2-zero"))
coldata <- coldata %>%
  filter(!is.na(condition))

# Extract B-values
data.meth <- dna.meth %>% 
  assay() %>%
  as.data.frame()%>%
  na.omit() 

# Select columns by metadata
selected_columns <- coldata$barcode
data.meth <- data.meth[, selected_columns]

# Purity data 
data.purity <- read.xlsx("TCGA_cpe_purity.xlsx", startRow = 3)
data.purity <- data.purity %>%
  select(-X8) %>%
  filter(Cancer.type %in% "BRCA")

# Select columns by purity
data.purity <- data.purity %>% 
  filter(Sample.ID %in% metadata$Sample.ID)

data.purity.high <- data.purity %>%
  filter(CPE > 0.6)

coldata <- coldata %>% 
  filter(sample %in% data.purity.high$Sample.ID)

selected_colums2 <- coldata$barcode

data.meth <- data.meth[, selected_colums2]


# Subset groups
colnames(data.meth) <- substr(colnames(data.meth), 1,16)
myCombat.HER2low <- data.meth %>%
  select(any_of(meta_HER2low$Sample.ID)) # 33 patients/ data HER2-low

myCombat.HER2zero <- data.meth %>%
  select(any_of(meta_HER2zero$Sample.ID)) # 15 patients/ data HER2-zero


metadata <- metadata %>%
  filter(Sample.ID %in% colnames(data.meth))

meta_HER2low <- metadata %>%
  filter(Group %in% "HER2-low")
meta_HER2zero <- metadata %>%
  filter(Group %in% "HER2-zero")


# Extract means and fold
HER2low.mean <- apply(myCombat.HERlow, 1, mean)
HER2zero.mean <- apply(dataTNBC, 1, mean)

fold <- HER2low.mean - HER2zero.mean 

# Compute statistical significance #
pvalue = NULL
tstat = NULL
for(i in 1 : nrow(myCombat.HER2low)) {
  x = myCombat.HER2low[i,]
  y = myCombat.HER2zero[i,]
  
  t = wilcox.test(as.numeric(x), as.numeric(y)) 
  pvalue[i] = t$p.value
  tstat[i] = t$statistic
}


all_probes <- as.data.frame(cbind(pvalue, HER2low.mean, HER2zero.mean, fold, abs(fold)))
all_probes <- all_probes %>%
  rownames_to_column(var = "IlmnID")
colnames(all_probes) <- c("IlmnID", "pvalue", "avg.HER2low", "avg.HER2zero", "fold", "abs(fold)")

# DMSs
all_probes <- all_probes %>%
  mutate(Meth_status = case_when(fold > 0.10 & pvalue < 0.05 ~ "Hypermethylated",
                                 fold < -0.10 & pvalue < 0.05 ~ "Hypomethylated",
                                 T ~ "NS")) %>%
  mutate(Meth_status = factor(Meth_status, levels = c("Hypomethylated", "NS", "Hypermethylated"))) 

table(all_probes$Meth_status)

probes_with_changes <- all_probes %>%
  filter(Meth_status %in% c("Hypermethylated", "Hypomethylated"))

probes_with_changes.hyper <- probes_with_changes %>%
  filter(Meth_status %in% "Hypermethylated")

probes_with_changes.hypo <- probes_with_changes %>%
  filter(Meth_status %in% "Hypomethylated")

myCombat.DMS <- data.meth[probes_with_changes$IlmnID, ]


# Loading manifest
library(readr)
manifest.450k <- read_tsv("C:/Users/idisb/Desktop/Datasets Epigenetics Cancer Lab/01 - Other Important information (Piedra roseta, Manifests, ...)/DNAmeth Manifests/HM450K.hg38.curated.txt")
manifest.450k <- manifest.450k %>%
  mutate(CGIposition = case_when(
    CGIposition %in% c("N_Shore", "S_Shore") ~ "Shore",
    CGIposition %in% c("S_Shelf", "N_Shelf") ~ "Shelf",
    TRUE ~ CGIposition
  ))


probes_and_manifest <- merge(all_probes, manifest.450k, by = 1)
probes_and_manifest_changes <- merge(probes_with_changes, manifest.450k, by = 1)

probes_and_manifest_changes.hyper <- probes_and_manifest_changes %>%
  filter(Meth_status %in% "Hypermethylated")
probes_and_manifest_changes.hypo <- probes_and_manifest_changes %>%
  filter(Meth_status %in% "Hypomethylated")



#### PLOTS AND ANALYSES ####

# Volcano Plot
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggrepel)

ggplot(all_probes, aes(x = fold, y = -log10(pvalue), color = Meth_status)) + 
  geom_point(size=1) + 
  xlim(-0.33,0.33) + ylim (0,6) + 
  theme_pubr() +
  scale_color_manual("Methylation status", values=c("limegreen", "grey", "red")) + 
  scale_x_continuous(breaks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3)) +
  geom_vline(xintercept = fold_cutoff, color = "black", linetype = "dashed", linewidth = 0.8) + 
  geom_vline(xintercept = -fold_cutoff, color = "black", linetype = "dashed", linewidth = 0.8) + 
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed", linewidth = 0.8) 
  


# Higher DNAm levels. Hypermethylated vs Hypomethylated

# Pvalue < 0.05 Yes/No
# fold > 0.1 or < -0.1
nrow(all_probes %>%
  filter(pvalue > 0.05 & fold > 0.1))
nrow(all_probes %>%
  filter(pvalue > 0.05 & fold < (-0.1)))

# Yes/Yes = 4357
#Yes/No = 930
#No/Yes = 3674
#Yes/No = 1477

a <- c(4357, 930)
b <- c(3674, 1477)
x <- rbind(a,b)
x <- as.table(x)
chisq.test(x)



# Enrichment analyses 
library(missMethyl)

# Gene Ontology Hypermethylated

set.seed(46)
go_meth.hyper <- gometh(probes_with_changes.hyper$IlmnID,
                  all.cpg = all_probes$IlmnID,
                  collection = c("GO"),
                  array.type = c("450K"),
                  plot.bias = FALSE,
                  prior.prob = TRUE,
                  anno = NULL,
                  equiv.cpg = FALSE,
                  fract.counts = TRUE,
                  genomic.features = c("TSS200", "TSS1500"),
                  sig.genes = TRUE)

go_meth.hyper.df <- go_meth.hyper %>%
  filter(FDR < 0.001) %>%
  filter(N > 3) %>%
  filter(ONTOLOGY %in% "BP")


go_meth.hyper.df <- go_meth.hyper.df %>%
  mutate(Ratio = (DE/N)*100) %>% # Ratio
  arrange((FDR))

top.10.hyper <- go_meth.hyper.df[1:10,]

top.10.hyper %>%
  mutate(TERM = fct_reorder(TERM, Ratio)) %>%
  ggplot(aes(x = TERM, y = Ratio, fill = FDR)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.6) +
  coord_flip() +
  xlab("") +
  theme_bw() + 
  scale_fill_gradient(low = "red", high = "blue")


# Biclusters
Bicluster1 <- read.xlsx("Biclusters/Bicluster 1 eMQTLs (Cell Cycle).xlsx")
Bicluster1 <- Bicluster1 %>%
  filter(Probe %in% all_probes$IlmnID)
Bicluster1_changes <- merge(Bicluster1, probes_with_changes, by = 1)

Bicluster2 <- read.xlsx("Biclusters/Bicluster 2 eMQTLs (Estrogen Response).xlsx")
Bicluster2 <- Bicluster2 %>%
  filter(Probe %in% all_probes$IlmnID)
Bicluster2_changes <- merge(Bicluster2, probes_with_changes, by = 1)

Bicluster3 <- read.xlsx("Biclusters/Bicluster 3 eMQTLs (EMT).xlsx")
Bicluster3 <- Bicluster3 %>%
  filter(Probe %in% all_probes$IlmnID)
Bicluster3_changes <- merge(Bicluster3, probes_with_changes, by = 1)

Bicluster4 <- read.xlsx("Biclusters/Bicluster 4 eMQTLs (Estrogen Response).xlsx")
Bicluster4 <- Bicluster4 %>%
  filter(Probe %in% all_probes$IlmnID)
Bicluster4_changes <- merge(Bicluster4, probes_with_changes, by = 1)

Bicluster5 <- read.xlsx("Biclusters/Bicluster 5 eMQTLs (Immune System).xlsx")
Bicluster5 <- Bicluster5 %>%
  filter(Probe %in% all_probes$IlmnID)
Bicluster5_changes <- merge(Bicluster5, probes_with_changes, by = 1)



## Hacemos el barplot %
B1.p <- (length(Bicluster1_changes$Probe)/length(Bicluster1$Probe)*100)
B2.p <- (length(Bicluster2_changes$Probe)/length(Bicluster2$Probe)*100)
B3.p <- (length(Bicluster3_changes$Probe)/length(Bicluster3$Probe)*100)
B4.p <- (length(Bicluster4_changes$Probe)/length(Bicluster4$Probe)*100)
B5.p <- (length(Bicluster5_changes$Probe)/length(Bicluster5$Probe)*100)


B.df <- data.frame(
  Bicluster = c("Bicluster 1 -  Cell cycle", "Bicluster 2 - Estrogen response", 
                "Bicluster 3 - EMT process", "Bicluster 4 - Estrogen response", "Bicluster 5 - Immune response"),
  Percentage = c(B1.p, B2.p, B3.p, B4.p, B5.p))

B.df$Bicluster <- factor(B.df$Bicluster, levels = rev(B.df$Bicluster))

library(ggpubr)
ggplot(B.df, aes(x = Percentage, y = Bicluster)) +
  geom_bar(stat = "identity", fill = "tomato", width = 0.8) +
  geom_text(aes(label = paste0(round(Percentage, 2), "%")),
            hjust = -0.2, size = 3, color = "black") +
  theme_pubr()



# Circus plot #
library(circlize)
DMR_hyper <- probes_and_manifest_changes.hyper %>%
  select(CpG_chrm, CpG_beg, CpG_end)

colnames(DMR_hyper) <- c("chr", "start", "end")

DMR_hypo <- probes_and_manifest_changes.hypo %>%
  select(CpG_chrm, CpG_beg, CpG_end)

colnames(DMR_hypo) <- c("chr", "start", "end")

circos.initializeWithIdeogram(chromosome.index = c(paste0("chr", 1:22), "chrX"))
circos.genomicDensity(DMR_hyper, col = c("red"), track.height = 0.1)
circos.genomicDensity(DMR_hypo, col = c("limegreen"), track.height = 0.1)


# To represent at Genome Browser
DMR_hyper.6 <- DMR_hyper %>%
  filter(chr == "chr6") %>%
  filter(start > 29602228 & end < 33409880) %>% # MHC region
  mutate(value = 0.5)

write_tsv(DMR_hyper.6, file = "BED HLA hypermethylation.bedGraph", col_names = FALSE)

circos.initializeWithIdeogram(chromosome.index = c(paste0("chr", 1:22), "chrX"))
circos.genomicDensity(DMR_hyper.6, col = c("red"), track.height = 0.1)