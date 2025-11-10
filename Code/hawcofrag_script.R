install.packages("crayon")
install.packages("readr")
install.packages("tidyverse")
install.packages("tibble")
install.packages("R.utils")
BiocManager::install("DESeq2")
library(apeglm)
library(tidyverse)
library(tibble)
library(readr)
library(crayon)
library(dplyr)
library(R.utils)
library(tidyverse)
library(DESeq2)

gunzip("~/Downloads/Fracy1_goinfo_FilteredModels2.tab.gz")
goinfo <- read_delim("~/Downloads/Fracy1_goinfo_FilteredModels2.tab", delim = "\t")
colnames(goinfo)[colnames(goinfo) == "#proteinId"] <- "proteinID"
colnames(goinfo)[colnames(goinfo) == "gotermId"] <- "GOtermID"
View(goinfo)
goinfo$proteinID <- as.character(goinfo$proteinID)
goinfo$GOtermID <- as.character(goinfo$GOtermID)

domaininfo <- read_delim("~/Downloads/fracy1_domaininfo_filteredmodels3.tab", delim = "\t")
colnames(domaininfo)[colnames(domaininfo) == "#proteinId"] <- "proteinID"
domaininfo$proteinID <- as.character(domaininfo$proteinID)
domaininfo$numHits <- as.character(domaininfo$numHits)
domaininfo$domainEnds <- as.character(domaininfo$domainEnds)
newannotations <- left_join(goinfo, domaininfo, by = "proteinID")
newannotations_merged <- left_join(merged, newannotations, by = "proteinID")

counts <- read_delim("~/Downloads/Lizzy/hawcofrag/merged_counts_sorted.txt", delim = "\t")
colnames(counts) <- c(
  "proteinID", "geneID",
  "sample1", "sample2", "sample3", "sample4", "sample5", "sample6",
  "sample7", "sample8", "sample9", "sample10", "sample11", "sample12"
)
kogannotations <- read_delim("~/Downloads/Lizzy/hawcofrag/Fracy1_koginfo_FilteredModels2.sorted1.txt", delim = "\t")
  kogannotations <- kogannotations[, -2]
  kogannotations <- kogannotations[-nrow(kogannotations),]
  colnames(kogannotations) <- c("proteinID", "kog_id", "kog_defline", "kog_class", "kog_group")
ecannotations <- read_delim("~/Downloads/Lizzy/hawcofrag/Fracy1_ecpathwayinfo_FilteredModels2.sorted.txt", delim = "\t")
  colnames(ecannotations) <- c("proteinID", "ec_number", "ec_definition", "ec_catalytic_activity", "ec_cofactors", "ec_associated_diseases", "ec_pathway", "ec_pathway_class", "ec_pathway_type")

ecannotations$proteinID <- as.character(ecannotations$proteinID)
kogannotations$proteinID <- as.character(kogannotations$proteinID)
counts$proteinID <- as.character(counts$proteinID)

merged <- left_join(counts, kogannotations, by = "proteinID")
merged <- left_join(merged, ecannotations, by = "proteinID")

merged_annotations <- merged[,-c(3:14)]
new_merged_annotations <- newannotations_merged[,-c(3:14)]

write.table(
  x = merged,
  file = "counts_annotated1.txt",
  sep = "\t", 
  row.names = FALSE 
)

counts_geneID <- counts[, -1] # counts with geneID only 
rows_to_remove <- c(3039, 13148 , 15887)  
counts_geneID <- counts_geneID[-rows_to_remove, ]
counts_geneID <- counts_geneID %>%
  remove_rownames() %>%
  column_to_rownames(var = names(.)[1])

samples <- colnames(counts_geneID) # Sample names from count matrix

coldata <- data.frame(
  row.names = samples,
  condition = c("noMn", "23Mn33Zn", "23Mn0.08Zn", "noMn", "23Mn33Zn", "noMn", "23Mn33Zn", "23Mn0.08Zn", "23Mn0.08Zn", "noMn33Zn", "noMn33Zn", "noMn33Zn")
)
coldata

dds <- DESeqDataSetFromMatrix(
  countData = counts_geneID,
  colData = coldata,
  design = ~ condition
)
dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

new_normalized_annotated <- inner_join(new_merged_annotations, normalized_counts_df, by = "geneID")
new_collapsed <- new_normalized_annotated %>%
  group_by(geneID) %>%
  summarise(across(everything(), ~ paste(unique(.x), collapse = ", ")))
write.csv(as.data.frame(new_collapsed), file = "new_counts_normalized_annotated_collapsed.csv", row.names = FALSE)


vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)  # this gives you the VST-normalized counts matrix
vst_normalized_counts <- as.data.frame(vsd_mat) %>%
  rownames_to_column(var = "geneID") 
new_merged_annotations$geneID <- as.character(new_merged_annotations$geneID)
vst_normalized_counts$geneID <- as.character(vst_normalized_counts$geneID)
vst_normalized_annotated <- inner_join(new_merged_annotations,
                                   vst_normalized_counts,
                                   by = "geneID")
vst_collapsed <- vst_normalized_annotated %>%
  group_by(geneID) %>%
  summarise(across(everything(), ~ paste(unique(.x), collapse = ", ")))

normalized_counts_df <- as.data.frame(normalized_counts) %>%
  rownames_to_column(var = "geneID") 
merged_annotations$geneID <- as.character(merged_annotations$geneID)
normalized_counts_df$geneID <- as.character(normalized_counts_df$geneID) 
normalized_annotated <- inner_join(merged_annotations,
                                   normalized_counts_df,
                                   by = "geneID")
collapsed <- normalized_annotated %>%
  group_by(geneID) %>%
  summarise(across(everything(), ~ paste(unique(.x), collapse = ", ")))
write.csv(as.data.frame(collapsed), file = "counts_normalized_annotated_collapsed.csv", row.names = FALSE)


collapsed_for_DE <- new_collapsed[,-c(28:39)]
write.csv(as.data.frame(collapsed_for_DE), file = "collapsed_for_DE_annotations.csv", row.names = FALSE)

# Visualization 
res1LFC <- lfcShrink(dds, coef="condition_23Mn33Zn_vs_23Mn0.08Zn", type="apeglm")
res_df1 <- as.data.frame(res1LFC)
res_df1$gene <- rownames(res_df1)
res_df1 <- res_df1[!is.na(res_df1$padj), ]
top_genes_res1 <- res_df1[order(res_df1$padj, -abs(res_df1$log2FoldChange)), ][1:10, ]
plotMA(res1LFC, alpha = 0.05, ylim=c(-13,13))
title("Differential Expression: high_Mn high_Zn vs high_Mn low_Zn")
 with(top_genes_res1
      , text(baseMean, log2FoldChange, labels = gene, pos = 3, cex = 0.8, col = "red"))
 list(top_genes_res1)

res1 <- results(dds, name ="condition_23Mn33Zn_vs_23Mn0.08Zn", alpha=0.05)
summary(res1)
sum(res1$padj < 0.05, na.rm=TRUE)
res1Ordered <- res1[order(res1$padj),]
write.csv(as.data.frame(res1Ordered), 
          file="highMnhighZn_highMnlowZn_DESeq.csv")
res1Ordered_df <- as.data.frame(res1Ordered) %>%
  rownames_to_column(var = "geneID")
res1_annotated <- inner_join(collapsed_for_DE, res1Ordered_df, by = "geneID")


res2LFC <- lfcShrink(dds, coef="condition_noMn_vs_23Mn0.08Zn", type="apeglm")
res_df2 <- as.data.frame(res2LFC)
res_df2$gene <- rownames(res_df2)
top_genes_res2 <- res_df2[order(res_df2$padj, -abs(res_df2$log2FoldChange)), ][1:10, ]
plotMA(res2LFC, alpha = 0.05, ylim=c(-13,13))
title("Differential Expression: no_Mn vs high_Mn low_Zn")
with(top_genes_res2
     , text(baseMean, log2FoldChange, labels = gene, pos = 3, cex = 0.8, col = "red"))
list(top_genes_res2)

res2 <- results(dds, name ="condition_noMn_vs_23Mn0.08Zn", alpha=0.05)
summary(res2)
sum(res2$padj < 0.05, na.rm=TRUE)
res2Ordered <- res2[order(res2$padj),]
write.csv(as.data.frame(res2Ordered), 
          file="noMn_highMnlowZn_DESeq.csv")
res2Ordered_df <- as.data.frame(res2Ordered) %>%
  rownames_to_column(var = "geneID")
res2_annotated <- inner_join(collapsed_for_DE, res2Ordered_df, by = "geneID")


res3LFC <- lfcShrink(dds, coef="condition_noMn33Zn_vs_23Mn0.08Zn", type="apeglm")
res_df3 <- as.data.frame(res3LFC)
res_df3$gene <- rownames(res_df3)
res_df3 <- res_df3[!is.na(res_df3$padj), ]
top_genes_res3 <- res_df3[order(res_df3$padj, -abs(res_df3$log2FoldChange)), ][1:10, ]
plotMA(res3LFC, alpha = 0.05, ylim=c(-13,13))
title("Differential Expression: no_Mn high_Zn vs high_Mn low_Zn")
with(top_genes_res3
     , text(baseMean, log2FoldChange, labels = gene, pos = 3, cex = 0.8, col = "red"))
list(top_genes_res3)

res3 <- results(dds, name ="condition_noMn33Zn_vs_23Mn0.08Zn", alpha=0.05)
summary(res3)
sum(res3$padj < 0.05, na.rm=TRUE)
res3Ordered <- res3[order(res3$padj),]
write.csv(as.data.frame(res3Ordered), 
          file="noMnhighZn_highMnlowZn_DESeq.csv")
res3Ordered_df <- as.data.frame(res3Ordered) %>%
  rownames_to_column(var = "geneID")
res3_annotated <- inner_join(collapsed_for_DE, res3Ordered_df, by = "geneID")


#getting rest of pairwise expressions
dds$condition <- relevel(dds$condition, ref = "noMn")
dds <- DESeq(dds)
resultsNames(dds)

res4LFC <- lfcShrink(dds, coef="condition_23Mn33Zn_vs_noMn", type="apeglm")
res_df4 <- as.data.frame(res4LFC)
res_df4$gene <- rownames(res_df4)
res_df4 <- res_df4[!is.na(res_df4$padj), ]
top_genes_res4 <- res_df4[order(res_df4$padj, -abs(res_df4$log2FoldChange)), ][1:10, ]
plotMA(res4LFC, alpha = 0.05, ylim=c(-13,13))
title("Differential Expression: high_Mn high_Zn vs no_Mn")
with(top_genes_res4
     , text(baseMean, log2FoldChange, labels = gene, pos = 3, cex = 0.8, col = "red"))
list(top_genes_res4)

res4 <- results(dds, name ="condition_23Mn33Zn_vs_noMn", alpha=0.05)
summary(res4)
sum(res4$padj < 0.05, na.rm=TRUE)
res4Ordered <- res4[order(res4$padj),]
write.csv(as.data.frame(res4Ordered), 
          file="highMnhighZn_noMn_DESeq.csv")
res4Ordered_df <- as.data.frame(res4Ordered) %>%
  rownames_to_column(var = "geneID")
res4_annotated <- inner_join(collapsed_for_DE, res4Ordered_df, by = "geneID")


res5LFC <- lfcShrink(dds, coef="condition_noMn33Zn_vs_noMn", type="apeglm")
res_df5 <- as.data.frame(res5LFC)
res_df5$gene <- rownames(res_df5)
res_df5 <- res_df5[!is.na(res_df5$padj), ]
top_genes_res5 <- res_df5[order(res_df5$padj, -abs(res_df5$log2FoldChange)), ][1:10, ]
plotMA(res5LFC, alpha = 0.05, ylim=c(-13,13))
title("Differential Expression: no_Mn high_Zn vs. no Mn")
with(top_genes_res5
     , text(baseMean, log2FoldChange, labels = gene, pos = 3, cex = 0.8, col = "red"))
list(top_genes_res5)

res5 <- results(dds, name ="condition_noMn33Zn_vs_noMn", alpha=0.05)
summary(res5)
sum(res5$padj < 0.05, na.rm=TRUE)
res5Ordered <- res5[order(res5$padj),]
write.csv(as.data.frame(res5Ordered), 
          file="noMnhighZn_noMn_DESeq.csv")
res5Ordered_df <- as.data.frame(res5Ordered) %>%
  rownames_to_column(var = "geneID")
res5_annotated <- inner_join(collapsed_for_DE, res5Ordered_df, by = "geneID")


#getting last pairwise expression 
dds$condition <- relevel(dds$condition, ref = "23Mn33Zn") 
dds <- DESeq(dds)
resultsNames(dds)

res6LFC <- lfcShrink(dds, coef="condition_noMn33Zn_vs_23Mn33Zn", type="apeglm")
res_df6 <- as.data.frame(res6LFC)
res_df6$gene <- rownames(res_df6)
res_df6 <- res_df6[!is.na(res_df6$padj), ]
top_genes_res6 <- res_df6[order(res_df6$padj, -abs(res_df6$log2FoldChange)), ][1:10, ]
plotMA(res6LFC, alpha = 0.05, ylim=c(-13,13))
title("Differential Expression: no_Mn high_Zn vs. high_Mn high_Zn") 
with(top_genes_res6
     , text(baseMean, log2FoldChange, labels = gene, pos = 3, cex = 0.8, col = "red"))
list(top_genes_res6)

res6 <- results(dds, name ="condition_noMn33Zn_vs_23Mn33Zn", alpha=0.05)
summary(res6)
sum(res6$padj < 0.05, na.rm=TRUE)
res6Ordered <- res6[order(res6$padj),]
write.csv(as.data.frame(res6Ordered), 
          file="noMnhighZn_highMnhighZn_DESeq.csv")
res6Ordered_df <- as.data.frame(res6Ordered) %>%
  rownames_to_column(var = "geneID")
res6_annotated <- inner_join(collapsed_for_DE, res6Ordered_df, by = "geneID")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("regionReport")
BiocManager::install("pcaExplorer")
library(pcaExplorer)
pcaExplorer(countmatrix = counts_geneID, coldata = coldata)

# converting EC number annotations to KO annotations 
BiocManager::install("KEGGREST")
library(KEGGREST)
get_annotation_for_ec <- function(ec) {
  res_link <- tryCatch({
    keggLink("ko", paste0("ec:", ec))
  }, error = function(e) NULL)
  
  if (is.null(res_link) || length(res_link) == 0) {
    return(data.frame(
      EC = ec,
      KO = NA,
      KO_name = NA,
      KO_definition = NA,
      pathways = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  ko_ids <- unique(unname(res_link))
  
  details_list <- lapply(ko_ids, function(ko) {
    kk <- keggGet(ko)
    entry <- kk[[1]]
    name <- if (!is.null(entry$NAME)) entry$NAME else NA
    definition <- if (!is.null(entry$DEFINITION)) entry$DEFINITION else NA
    pathways <- if (!is.null(entry$PATHWAY)) {
      paste(names(entry$PATHWAY), entry$PATHWAY, collapse = "; ")
    } else {
      NA
    }
    data.frame(
      EC = ec,
      KO = ko,
      KO_name = name,
      KO_definition = definition,
      pathways = pathways,
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, details_list)
}
# --- Progress bar wrapper ---
ec_vec <- gsub("^ec:", "", trimws(new_collapsed$ec_number))
ec_vec <- ec_vec[!is.na(ec_vec) & ec_vec != ""]

pb <- txtProgressBar(min = 0, max = length(ec_vec), style = 3)

results_list <- vector("list", length(ec_vec))
for (i in seq_along(ec_vec)) {
  results_list[[i]] <- get_annotation_for_ec(ec_vec[i])
  setTxtProgressBar(pb, i)
  Sys.sleep(0.2)  # optional pause to be gentle to KEGG
}
close(pb)

# adding KO annotations to all annotation tables
ko_annotation_table <- do.call(rbind, results_list)
colnames(ko_annotation_table) <- c("ec_number", "ko_ID", "ko_name", "ko_definition", "ko_pathway")
ko_annotation_table$ec_number <- as.character(ko_annotation_table$ec_number)
ko_annotation_table$ec_number <- trimws(ko_annotation_table$ec_number)
ko_annotation_table$ec_number[ko_annotation_table$ec_number %in% c("", "NA", "n/a", "N/A")] <- NA
ko_annotation_table <- ko_annotation_table[!is.na(ko_annotation_table$ec_number), ]
  ko_annotation_table_collapsed <- ko_annotation_table %>%
  group_by(ec_number) %>%
  summarise(across(everything(), ~ paste(unique(.x), collapse = ", ")))

new_collapsed_ko <- left_join(new_collapsed, ko_annotation_table_collapsed, by = "ec_number")
new_collapsed_ko <- new_collapsed_ko %>%
  select(
    geneID,
    proteinID,
    ec_number,
    ko_ID, ko_name, ko_definition, ko_pathway, 
    everything()  
  )
write.csv(as.data.frame(new_collapsed_ko), file="FULLannotations_fcyl_normalized.csv")

vst_collapsed_ko <- left_join(vst_collapsed, ko_annotation_table_collapsed, by = "ec_number")
vst_collapsed_ko <- vst_collapsed_ko %>%
  select(
    geneID,
    proteinID,
    ec_number,
    ko_ID, ko_name, ko_definition, ko_pathway,
    everything() 
  )
write.csv(as.data.frame(vst_collapsed_ko), file="VST_FULLannotations_fcyl.csv")

res1_annotated_ko <- left_join(res1_annotated, ko_annotation_table_collapsed, by = "ec_number")
res1_annotated_ko <- res1_annotated_ko %>%
  select(
    geneID,
    proteinID,
    ec_number,
    ko_ID, ko_name, ko_definition, ko_pathway,
    everything()
  )
write.csv(as.data.frame(res1_annotated_ko), file="FULLannotations_fcyl_res1DE.csv")

res2_annotated_ko <- left_join(res2_annotated, ko_annotation_table_collapsed, by = "ec_number")
res2_annotated_ko <- res2_annotated_ko %>%
  select(
    geneID,
    proteinID,
    ec_number,
    ko_ID, ko_name, ko_definition, ko_pathway,
    everything() 
  )
write.csv(as.data.frame(res2_annotated_ko), file="FULLannotations_fcyl_res2DE.csv")

res3_annotated_ko <- left_join(res3_annotated, ko_annotation_table_collapsed, by = "ec_number")
res3_annotated_ko <- res3_annotated_ko %>%
  select(
    geneID,
    proteinID,
    ec_number,
    ko_ID, ko_name, ko_definition, ko_pathway, 
    everything() 
  )
write.csv(as.data.frame(res3_annotated_ko), file="FULLannotations_fcyl_res3DE.csv")

res4_annotated_ko <- left_join(res4_annotated, ko_annotation_table_collapsed, by = "ec_number")
res4_annotated_ko <- res4_annotated_ko %>%
  select(
    geneID,
    proteinID,
    ec_number,
    ko_ID, ko_name, ko_definition, ko_pathway, 
    everything() 
  )
write.csv(as.data.frame(res4_annotated_ko), file="FULLannotations_fcyl_res4DE.csv")

res5_annotated_ko <- left_join(res5_annotated, ko_annotation_table_collapsed, by = "ec_number")
res5_annotated_ko <- res5_annotated_ko %>%
  select(
    geneID,
    proteinID,
    ec_number,
    ko_ID, ko_name, ko_definition, ko_pathway, 
    everything() 
  )
write.csv(as.data.frame(res5_annotated_ko), file="FULLannotations_fcyl_res5DE.csv")

res6_annotated_ko <- left_join(res6_annotated, ko_annotation_table_collapsed, by = "ec_number")
res6_annotated_ko <- res6_annotated_ko %>%
  select(
    geneID,
    proteinID,
    ec_number,
    ko_ID, ko_name, ko_definition, ko_pathway, 
    everything() 
  )
write.csv(as.data.frame(res6_annotated_ko), file="FULLannotations_fcyl_res6DE.csv")
