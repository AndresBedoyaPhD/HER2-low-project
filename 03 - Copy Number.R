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


query_copy_number <- GDCquery(project = "TCGA-BRCA",
                              data.category = "Copy Number Variation",
                              access = "open",
                              data.type = "Gene Level Copy Number",
                              barcode = metadata$Sample.ID,
                              sample.type = "Primary Tumor")

data_copy <- getResults(query_copy_number)

# Download data 
GDCdownload(query_copy_number)

copy.number <- GDCprepare(query_copy_number) 

data.copy <- copy.number %>% 
  assay() %>%
  as.data.frame() %>%
  na.omit()


# Categorical classification
# 0 = Loss
# 1 = Heterocigozity
# 2 = Normal
# > 2  = Gain

metadata_filtered <- metadata %>%
  filter(Sample.ID %in% colnames(data.copy))%>%
  select(Patient.ID, Sample.ID, Age, CPE, Histology_summary, BC_Subtype_HER2) %>%
  mutate(Condition = case_when(BC_Subtype_HER2 %in% "4.1_TNBC_HER2_0" ~ "HER2-zero",
                               BC_Subtype_HER2 %in% "4.2_TNBC_HER2low" ~ "HER2-low")) %>%
  select(-BC_Subtype_HER2)

colnames(metadata_filtered) <- c("Patient.ID", "Sample.ID", "Age", "CPE", "Cancer.Type", "Condition")

gene_metadata <- as.data.frame(rowData(copy.number))

data.copy2 <- data.copy %>%
  rownames_to_column(var = "gene_id") %>%
  inner_join(gene_metadata, by = "gene_id") %>%
  select(-gene_id) %>%
  select(gene_name, everything())

HER2.copy <- data.copy2 %>%
  filter(gene_name %in% "ERBB2") %>%
  pivot_longer(-gene_name,names_to = "Sample.ID", values_to = "ERBB2_copy")%>%
  inner_join(metadata_filtered, by="Sample.ID")

HER2.copy <- HER2.copy %>%
  mutate(ERBB2_category = case_when(
    ERBB2_copy < 2 ~ "Loss",
    ERBB2_copy == 2 ~ "Diploid",
    ERBB2_copy > 2 ~ "Gain"
    ))

HER2.copy$ERBB2_category <- factor(HER2.copy$ERBB2_category, levels = c("Loss","Diploid", "Gain"))
round(prop.table(table(HER2.copy$ERBB2_category, HER2.copy$Condition), margin = 2), 2) * 100


HER2.copy2 <- HER2.copy %>%
  mutate(value = case_when(Condition == "HER2-zero" & ERBB2_category == "Loss" ~ 35,
                           Condition == "HER2-zero" & ERBB2_category == "Diploid" ~ 27,
                           Condition == "HER2-zero" & ERBB2_category == "Gain" ~ 38,
                           Condition == "HER2-low" & ERBB2_category == "Loss" ~ 14,
                           Condition == "HER2-low" & ERBB2_category == "Diploid" ~ 44,
                           Condition == "HER2-low" & ERBB2_category == "Gain" ~ 42))

HER2.copy2 <- HER2.copy2 %>%
  select(Condition, ERBB2_category, value) %>%
  filter(!duplicated(paste(Condition, ERBB2_category, value))) %>%
  arrange(ERBB2_category)


HER2.copy2$Condition <- factor(HER2.copy2$Condition, levels = c("HER2-low", "HER2-zero"))

ggplot(HER2.copy2, aes(y=ERBB2_category, x = value, fill = Condition))+ 
  geom_bar(position = "dodge", stat = "identity", width = 0.6) +
  theme_pubr() + 
  scale_fill_manual(values = c("HER2-low" = "#7CEBB0", "HER2-zero" = "#884B99")) +
  geom_text(aes(label = round(value, 2)),
            position = position_dodge(width = 0.6), 
            vjust = 0.5, hjust = -0.5, color = "black") +
  labs(x = "Percentage of ERBB2 copy number for each condition", y = "ERBB2 copy category")
  

chisq.test(table(HER2.copy$ERBB2_category, HER2.copy$Condition))
