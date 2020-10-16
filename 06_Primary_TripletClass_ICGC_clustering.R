
## clustering ICGC samples with
## Triplet genes from my samples

setwd("~/Desktop/ServerView/MYC/scripts/MycHdaci_code/scripts")

library(tidyverse)
library(dendextend)

library(ComplexHeatmap)
library(circlize)

all_gene_classes_exp <- readRDS("../data/Primary_AllClasses_matched_pid_exp.rds")

icgc_tpm <- data.table::fread("~/Desktop/ServerView/MYC/results/data_290917/rnaseq/kallisto/kallisto_downstream/ICGC_MB_RNAseq_normalized_tpm_matrix_subgroupnames.txt", header = TRUE)

## required PID with MYC exp and status

# ICGC_MB131 515.22459 - amp
# ICGC_MB291 479.12427 - amp
# ICGC_MB89 325.01230 - noamp
# ICGC_MB260 152.66733 - amp
# ICGC_MB111  86.64376 - noamp

## MYC NON OE
# ICGC_MB146 1.9030352 - noamp
# ICGC_MB166 1.5639215 - noamp
# ICGC_MB165 1.0837035 - noamp
# ICGC_MB168 0.4280805 - noamp
# ICGC_MB134 0.4010341 - noamp

sample_status <- data.frame(PID = c("ICGC_MB131", "ICGC_MB291", "ICGC_MB89", "ICGC_MB260", "ICGC_MB111",
                                    "ICGC_MB146", "ICGC_MB166", "ICGC_MB165", "ICGC_MB168", "ICGC_MB134"),
                            AMP = c(1, 1, 0, 1, 0, 0, 0, 0, 0, 0),
                            EXP = c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0),
                            AMP_CLR = c("#FFC000FF", "#FFC000FF", "lightgray", "#FFC000FF", "lightgray",
                                        "lightgray", "lightgray", "lightgray", "lightgray", "lightgray"),
                            EXP_CLR = c(rep("#FFC000FF", 5), c(rep("lightgray", 5))))

req_icgc_tpm <- icgc_tpm %>% 
  dplyr::select(matches("target|MB131|MB291|MB89|MB260|MB111|MB146|MB166|MB165|MB168|MB134"))

## get triplet genes from ICGC PID
req_icgc_tpm_triplet <- all_gene_classes_exp %>% 
  dplyr::filter(CLASS == "MYC_HDAC2_H3K27ac") %>% 
  dplyr::filter(MEAN > 0.5) %>% 
  dplyr::select(target_id) %>% 
  merge(req_icgc_tpm, by = "target_id")

## triplet bound & expressed from my samples
colnames(req_icgc_tpm_triplet) <- gsub("G3_", "", colnames(req_icgc_tpm_triplet))
req_icgc_tpm_triplet <- req_icgc_tpm_triplet %>% 
  tibble::column_to_rownames(var = "target_id")

req_icgc_tpm_triplet_lg <- log2(req_icgc_tpm_triplet+1)

## top 500 - all genes
gene_num <- c(200, 500, 1000, 1500, 2000, 3000, 4000)

pdf("../figures/Primary_TripletExpressed_G3ICGC_clustering.pdf", height = 4.5, width = 7.5)

for(i in 1:length(gene_num)){
  
  my_gene_num <- gene_num[i]
  
  rv <- genefilter::rowVars(req_icgc_tpm_triplet_lg)
  idx <- order(-rv)[1:my_gene_num]
  idx_mat <- req_icgc_tpm_triplet_lg[idx, ]
  
  my_dend <- idx_mat %>% t() %>% 
    dist %>% hclust %>% 
    as.dendrogram() 
  
  ## labels
  lab_order <- labels(my_dend)
  
  my_clr_order <- sample_status[match(lab_order, sample_status$PID), ]
  
  ## bars
  myc_amp_bar <- my_clr_order$AMP_CLR %>% as.character()
  myc_exp_bar <- my_clr_order$EXP_CLR %>% as.character()
  
  the_bars <- cbind(myc_amp_bar, myc_exp_bar)
  
  par(mar = c(10,8,1,1))
  plot(my_dend, cex.lab = 10)
  title(paste0("Triplet expressed, n = ", my_gene_num))
  colored_bars(colors = the_bars, dend = my_dend, 
               rowLabels = c("MYC Amplification", "MYC Expression"),
               sort_by_labels_order = FALSE, cex.lab = 26)
  
}

dev.off()
  
## label
clr_bar <- data.frame(Status = c("True", "False"),
           Binary = c(1, 0))

pdf("../figures/Primary_TripletExpressed_G3ICGC_clustering_colorLegend.pdf", height = 1, width = 1.8)

clr_bar %>% 
  reshape2::melt() %>% 
  ggplot(aes(variable, Status, fill = factor(value)))+
  geom_tile()+
  scale_fill_manual(values = c("0" = "lightgray", "1" = "#FFC000FF"))+
  scale_y_discrete(position = "right")+
  theme_classic(base_size = 24)+
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black", face = "bold"))

dev.off()


## complexHeatmap version

hc <- idx_mat %>% t() %>% 
  dist %>% hclust

anno_df <- hc$labels %>% as.data.frame() %>% 
  dplyr::rename(PID = ".")

## amp  no-amp
## 291  134
## 131  166
## 89   146
## 111  165
## 260  168

#plot(hc)
#sample_status[match(hc$labels, sample_status$PID), ]

my_order <- c("ICGC_MB291", "ICGC_MB131","ICGC_MB89","ICGC_MB111","ICGC_MB260",
              "ICGC_MB134","ICGC_MB166","ICGC_MB146","ICGC_MB165","ICGC_MB168")

sample_status2 <- sample_status
sample_status2$amp_binary <- ifelse(sample_status2$AMP == 1, "yes", "no")
sample_status2$exp_binary <- ifelse(sample_status2$EXP == 1, "yes", "no")

amp_clr <- sample_status[match(my_order, sample_status$PID), ] %>% 
  dplyr::pull(AMP_CLR) %>% as.character()
names(amp_clr) <- my_order

exp_clr <- sample_status[match(my_order, sample_status$PID), ] %>% 
  dplyr::pull(EXP_CLR) %>% as.character()
names(exp_clr) <- my_order


amp_anno <- HeatmapAnnotation(`MYC expression` = anno_df$PID,
                              `MYC amplification` = anno_df$PID,
                  col=list(`MYC amplification` = amp_clr,
                           `MYC expression` = exp_clr),
                  annotation_height = unit(c(5, 5), "cm"),
                  annotation_name_gp = gpar(fontsize = 16, face = "bold"),
                  border = TRUE,
                  show_legend = FALSE,
                  gap = unit(c(0.15), "cm"))

tt <- matrix(nc = 10, nr = 0)
colnames(tt) <- hc$labels

pdf('../figures/Primary_TripletExpressed_G3ICGC_clustering.pdf', height = 6, width = 7.5)

Heatmap(tt, cluster_columns = hc,
        column_split = 2,
        column_dend_height = unit(7, "cm"),
        top_annotation = amp_anno,
        column_title = "Triplet expressed genes (n=1000)",
        column_names_side = "top")

dev.off()

##
# ICGC_MB131 515.22459 - amp
# ICGC_MB291 479.12427 - amp
# ICGC_MB89 325.01230 - noamp
# ICGC_MB260 152.66733 - amp
# ICGC_MB111  86.64376 - noamp

## MYC NON OE
# ICGC_MB146 1.9030352 - noamp
# ICGC_MB166 1.5639215 - noamp
# ICGC_MB165 1.0837035 - noamp
# ICGC_MB168 0.4280805 - noamp
# ICGC_MB134 0.4010341 - noamp

anno_df2 <- anno_df
anno_df2$myc_exp <- c(86.64376, 515.22459, 0.4010341, 1.9030352, 1.0837035, 1.5639215,
                      0.4280805, 152.66733, 479.12427, 325.01230)

tt_pts <- HeatmapAnnotation(#`MYC expression` = anno_df2$PID,
                            `MYC amplification` = anno_df2$PID,
                            `MYC expression (TPM)` = anno_points(anno_df2$myc_exp, size = unit(4, "mm")),
                            col=list(`MYC amplification` = amp_clr),
                            annotation_height = unit(c(0.7, 4), "cm"),
                            annotation_name_gp = gpar(fontsize = 16, face = "bold"),
                            border = TRUE,
                            show_legend = FALSE,
                            gap = unit(c(0.15), "cm"))


pdf('../figures/Primary_TripletExpressed_G3ICGC_clustering_exppoints.pdf', height = 7, width = 9)

Heatmap(tt, cluster_columns = hc,
        column_split = 2,
        column_dend_height = unit(7, "cm"),
        top_annotation = tt_pts,
        column_title = "Triplet expressed genes (n=1000)",
        column_names_side = "top")

dev.off()

##
load("../Myc_Hdac2_allObjects_simpleTheme2.RData")

