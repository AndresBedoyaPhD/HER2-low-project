library(openxlsx)
library(tidyverse)
library(TCGAbiolinks)
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

#### SINGLE NUCLEOTIDE VARIATIONS #### 
query_mutations.HER2low <- GDCquery(project = "TCGA-BRCA",
                                    data.category = "Simple Nucleotide Variation",
                                    access = "open",
                                    data.type = "Masked Somatic Mutation",
                                    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
                                    barcode = meta_HER2low$Sample.ID)


query_mutations.HER2zero <- GDCquery(project = "TCGA-BRCA",
                                 data.category = "Simple Nucleotide Variation",
                                 access = "open",
                                 data.type = "Masked Somatic Mutation",
                                 workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
                                 barcode = meta_HER2zero$Sample.ID)

# Download data 
GDCdownload(query_mutations.HER2low)
GDCdownload(query_mutations.HER2zero)

# Processing data
maf.HER2low <- GDCprepare(query_mutations.HER2low)
maf.HER2zero <- GDCprepare(query_mutations.HER2zero)

# We filter to match the patients we have in metadata
maf.HER2low$Tumor_Sample_Barcode <- substr(maf.HER2low$Tumor_Sample_Barcode, 1, 16)
maf.HER2zero$Tumor_Sample_Barcode <- substr(maf.HER2zero$Tumor_Sample_Barcode, 1, 16)

maf.HER2low.filt <- maf.HER2low %>%
  filter(Tumor_Sample_Barcode %in% meta_HER2low$Sample.ID)
maf.HER2zero.filt <- maf.HER2zero %>%
  filter(Tumor_Sample_Barcode %in% meta_HER2zero$Sample.ID)


# MAF object and analyses
library(maftools)
library(DT)

maf.HER2low.final <- read.maf(maf = maf.HER2low.filt)
maf.HER2zero.final <- read.maf(maf = maf.HER2zero.filt)

# Match medata to know the patients that we are comparing
meta.HER2low.filtered <- metadata %>%
  filter(Sample.ID %in% maf.HER2low.filt$Tumor_Sample_Barcode)

meta.HER2zero.filtered <- metadata %>%
  filter(Sample.ID %in% maf.HER2zero.filt$Tumor_Sample_Barcode)

# ANALYSIS #
# Summary
plotmafSummary(maf = maf.HER2low.final, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE) # HER2low
plotmafSummary(maf = maf.HER2zero.final, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE) # TNBC

# Oncoplots
oncoplot(maf = maf.HER2low.final, top = 10, removeNonMutated = TRUE)
oncoplot(maf = maf.HER2zero.final, top = 10, removeNonMutated = TRUE)


coOncoplot(maf.HER2zero.final, maf.HER2low.final,
           m1Name = "HER2-zero" , m2Name =  "HER2-Low")

coOncoplot(maf.HER2zero.final,maf.HER2low.final,
           m1Name = "HER2-zero" , m2Name = "HER2-Low",
           genes = c("CARM1", "EGFR", "ERBB2", "ERBB3", "ERBB4", "ESR1", "GRIP1", "HRAS", "IL6", "IL6R", "IL6ST", "MAP2K1", "MAPK1", "MAPK3", "PIK3CG", "PIK3R1", "RAF1", "SHC1", "SOS1", "STAT3"), 
           removeNonMutated = FALSE, keepGeneOrder = TRUE)


# Tumor mutational burden
tmb.HER2zero <- tmb(maf.HER2zero.final) %>%
  mutate(Condition = "HER2-zero")
tmb.HER2low <- tmb(maf.HER2low.final) %>%
  mutate(Condition = "HER2-low")

# Assess normality
shapiro.test(tmb.HER2low$total_perMB)
shapiro.test(tmb.HER2zero$total_perMB)

# Wilcoxon test
wilcox.test(tmb.HER2low$total_perMB, tmb.HER2zero$total_perMB)

tmb_all <- rbind(tmb.HER2low, tmb.HER2zero) %>%
  mutate(Condition = factor(Condition, levels = c("HER2-zero", "HER2-low")))

library(ggpubr)
ggplot(tmb_all, aes(x = Condition, y = total_perMB_log, fill = Condition)) +
  geom_violin(width = 0.6) +
  geom_boxplot(width = 0.3, outlier.shape = NA, fill = "white") +
  geom_jitter(color="black", size=1, width = 0.1) + 
  scale_fill_manual(values = c("HER2-zero" = "#884B99", "HER2-low" = "#7CEBB0")) +
  theme_pubr() +
  theme(legend.position = "none")
