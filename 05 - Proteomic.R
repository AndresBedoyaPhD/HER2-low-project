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


metadata$Color <- NA
metadata$Color[which(metadata$Group == "HER2low")] = "orange"
metadata$Color[which(metadata$Group == "TNBC")] = "black"


# query proteomic data
query_proteomics <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Proteome Profiling",
  data.type = "Protein Expression Quantification",
  access = "open",
  barcode = metadata$Sample.ID)

query_proteomics_results <- getResults(query_proteomics)

#Normalization Steps:
  # Calculate the median for each protein across all samples: Compute the median signal intensity for each protein across all samples to establish a baseline signal level for that protein across experimental conditions.

  # Subtract the median from each protein's values: Subtract the median obtained in step 1 from the signal intensities for that protein in each sample. This corrects for signal intensity differences across samples due to technical or non-specific biological factors.

  # Calculate the median for each sample across all proteins: Compute the median of signal intensities across all proteins within each sample, providing a baseline signal level for each sample across all analyzed proteins.

  # Subtract the median from each sample's values: Subtract the median from step 3 from the signal intensities of each protein in the respective sample. This adjusts for signal intensity variations among proteins within a specific sample due to technical or non-specific biological factors.
# more info: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/RPPA_intro/


GDCdownload(query_proteomics)
proteins <- GDCprepare(query_proteomics, summarizedExperiment = T)

# Annotations: https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files
protein.metadata <- read_tsv("TCGA_antibodies_descriptions.gencode.v36.tsv")

# Match
protein.tot <- inner_join(protein.metadata, proteins, by = "peptide_target")

protein.tot2 <- protein.tot %>%
  select(peptide_target, gene_name, contains("TCGA")) %>%
  separate_rows(gene_name) %>%
  na.omit()


# Filtering metadata by samples we have in matrix
metadata <- metadata %>%
  filter(Sample.ID %in% colnames(proteins))



# HER2 proteomics
HER2 <- protein.tot2 %>%
  filter(gene_name %in% "ERBB2") %>%
  select(-gene_name) %>%
  pivot_longer(-peptide_target, names_to = "Sample.ID", values_to = "Protein") %>%
  inner_join(metadata, by = "Sample.ID") %>%
  mutate(Group = factor(Group, levels = c("HER2-zero", "HER2-low")))


ggplot(HER2, aes(x = Group, y = Protein, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = c("HER2-low" = "#7CEBB0", "HER2-zero" = "#884B99")) +
  geom_jitter(color="black", size=1, width = 0.1) +
  facet_wrap(~peptide_target) +
  theme_pubr()


# statistics
prot.HER2low <- protein.tot2 %>%
  filter(gene_name %in% "ERBB2") %>%
  select(-gene_name) %>%
  filter(!duplicated(peptide_target)) %>%
  column_to_rownames(var = "peptide_target") %>%
  select(any_of(meta_HER2low$Sample.ID)) %>% # 29 patients/ data TNBC
  t() %>%
  as.data.frame()


prot.HER2zero <- protein.tot2 %>%
  filter(gene_name %in% "ERBB2") %>%
  select(-gene_name) %>%
  filter(!duplicated(peptide_target)) %>%
  column_to_rownames(var = "peptide_target") %>%
  select(any_of(meta_HER2zero$Sample.ID)) %>% # 29 patients/ data TNBC
  t() %>%
  as.data.frame()

# It's normal distribution
shapiro.test(prot.HER2low$HER2)
shapiro.test(prot.HER2zero$HER2)

# Not normal distribution
shapiro.test(prot.HER2low$HER2_pY1248) 
shapiro.test(prot.HER2zero$HER2_pY1248)

# HER2
t.test(prot.HER2low$HER2, prot.HER2zero$HER2)

# pHER2 Y1248
wilcox.test(prot.HER2low$HER2_pY1248, prot.HER2zero$HER2_pY1248