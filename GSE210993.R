# RNA seq

setwd("E:/MY PROJECT")
dir()

library(tidyverse)
library(ggplot2)

# Data Reprocessing

GSE210993 <- read.delim("GSE210993_IM.expression_TPM_matrix.103.txt", sep = "\t")
ncol(GSE210993)
colnames(GSE210993)

table(GSE210993$Gene_Type)
annot_tabel <- GSE210993[,c(1,7,8)]
view(annot_tabel)
write.csv(annot_tabel, "Annotation Tabel.csv")

samples <- colnames(GSE210993)
class(samples)
table(samples)

##########
# Metadata Preparation

group <- rep(NA, length(samples))
group[grepl("IM", samples)] <- "IM"
group[grepl("AO|BO|AgO|AlO|BgO|BlO|BfO|BfF|AF|BgF|BlF|AlF|AgF|BF", samples)] <- "Gastric"
group[grepl("NileO|NCaeO|NAsCO|NileF", samples)] <- "Intestinal"

metadata <- data.frame(
  Samples = samples,
  Group = factor(group)
)
class(metadata)

metadata <- metadata[c(-1:-13), ]
colnames(metadata) <- c("Sample_Name", "Sample_Type")
view(metadata)
write.csv(metadata, "MetaData.csv")


table(metadata$Sample_Type, useNA = "always")
view(metadata)
colnames(metadata)

GSE210993_subset <- GSE210993 |>
  select(any_of(metadata$Sample_Name), Gene_ID) |>
  column_to_rownames(var = "Gene_ID")

length(GSE210993$Gene_Name) == nrow(GSE210993_subset)
all(metadata$Sample_Name == colnames(GSE210993_subset))
identical(metadata$Sample_Name,
          colnames(GSE210993_subset))

view(GSE210993_subset)
ncol(GSE210993_subset)
nrow(GSE210993_subset)
dim(GSE210993_subset)

write.csv(GSE210993_subset, "Data Subset.csv")

#########
# Box plot

GSE210993_subset_long <- GSE210993_subset |>
  rownames_to_column("Gene_ID") |>
  pivot_longer(cols = -Gene_ID, names_to = "Sample", values_to = "Expression")
view(GSE210993_subset_long)

B1 <- ggplot(GSE210993_subset_long, aes(Sample, Expression)) +
    geom_boxplot(fill = "steelblue", alpha = 0.7) +
    theme_minimal()

pdf("Box_Plot_Subset.pdf", width = 8, height = 6)
B2 <- boxplot(
    log2(GSE210993_subset + 1),
    outline = FALSE,
    las = 2,
    main = "Boxplot Of log2(TPM + 1)"
  )
dev.off()

hist(GSE210993_subset[,1], breaks = 1000, xlim = c(0,50))
hist(GSE210993_subset[,2], breaks = 1000, xlim = c(0,50))
hist(GSE210993_subset[,3], breaks = 1000, xlim = c(0,50))

##########
# Gene Filtering

GSE210993_subset_filter_1 <- GSE210993_subset[which(rowSums(GSE210993_subset) != 0),]
nrow(GSE210993_subset_filter_1)


GSE210993_subset_filter_2 <- GSE210993_subset_filter_1[which(rowSums(GSE210993_subset_filter_1) > 10),]
nrow(GSE210993_subset_filter_2)

write.csv(GSE210993_subset_filter_2, "Data Subset Filtered.csv")

########
# Subset Gene Types

view(GSE210993_subset_filter_2)
view(annot_tabel)

annot_code <- GSE210993_subset_filter_2 |>
  rownames_to_column(var = "Gene_ID") |>
  inner_join(annot_tabel, by = "Gene_ID")

view(annot_code)
colnames(annot_code)
nrow(annot_code)
identical(
  rownames(GSE210993_subset_filter_2),
  annot_code[,1]
)

table(annot_code$Gene_Type)

ProteinCoding <- subset(annot_code, annot_code$Gene_Type == "protein_coding")
Pseudogene <- subset(annot_code, annot_code$Gene_Type == "pseudogene")
lincRNA <- subset(annot_code, annot_code$Gene_Type == "lincRNA")
antisense <- subset(annot_code, annot_code$Gene_Type == "antisense")

write.csv(ProteinCoding, "protein_coding.csv")
write.csv(Pseudogene, "pseudogene.csv")
write.csv(lincRNA, "lincRNA.csv")
write.csv(antisense, "antisense.csv")


#######
# Differential Expression Analysis using limma

view(GSE210993_subset_filter_2)
colnames(GSE210993_subset_filter_2)

# Step 1: Convert To Matrix  

step_1 <- as.matrix(GSE210993_subset_filter_2)
mode(step_1) <- "numeric"

# Step 2: Logarithm
step_2 <- log2(step_1 + 1)

# Step 3: Make Group
step_3 <- factor(metadata$Sample_Type)
table(step_3)

# step 4: Limma Library
library(limma)

# Step 5: Design Matrix
design <- model.matrix( ~0 + step_3)
colnames(design) <- levels(step_3)
design

# Step 6: Model
fit <- lmFit(step_2, design)


# step 7: Explain Contrast
contrast.matrix <- makeContrasts(
  IMvsGastric = 
    IM - Gastric,
  IMvsIntestinal =
    IM - Intestinal,
  GastricvsIntestinal =
    Gastric - Intestinal,
levels = design
)

# Step 8: Execute Contrast
fit2 <- contrasts.fit(
  fit,
  contrast.matrix
)

# Step 9: Empirical Bayes
fit2 <- eBayes(fit2)

# Step 10: Extract GEG

# Im vs Gastric
deg_Im_Gastric <- topTable(
  fit2,
  coef = "IMvsGastric",
  number = Inf
)

write.csv(deg_Im_Gastric, "DEG_IM_Gastric.csv")

# Im vs Intestinal
deg_Im_Intestinal <- topTable(
  fit2,
  coef = "IMvsIntestinal",
  number = Inf
)

write.csv(deg_Im_Intestinal, "DEG_IM_Intestinal.csv")

# Gastric vs Intestinal
deg_Gastric_Intestinal <- topTable(
  fit2,
  coef = "GastricvsIntestinal",
  number = Inf
)

write.csv(deg_Gastric_Intestinal, "DEG_Gastric_Intestinal.csv")

# Step 11:Filter Gene

deg_Im_Gastric_filtered <- subset(
  deg_Im_Gastric,
  adj.P.Val < 0.05 &
  abs(logFC) > 1
)

write.csv(deg_Im_Gastric_filtered, "DEG_IM_Gastric_Filtered.csv")

deg_Im_Intestinal_filtered <- subset(
  deg_Im_Intestinal,
  adj.P.Val < 0.05 &
  abs(logFC) > 1
)

write.csv(deg_Im_Intestinal_filtered, "DEG_IM_Intestinal_Filtered.csv")

deg_Gastric_Intestinal_filtered <- subset(
  deg_Gastric_Intestinal,
  adj.P.Val < 0.05 &
    abs(logFC) > 1
)

write.csv(deg_Gastric_Intestinal_filtered, "DEG_Gastric_Intestinal_Filtered.csv")

