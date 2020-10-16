## source: /b06x-lsdf/infra5-lsdf/mbHDAC2/scripts/data_290917/Figures/FeatureDistribution/featureDistribution_plot.R

## Feature distribution plot

setwd("~/Desktop/ServerView/MYC/scripts/MycHdaci_code/scripts")

library(dplyr)
library(ChIPseeker)
library(ggpubr)

library(ReactomePA)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)

library(cowplot)
library(ggplot2)
library(gtable)
library(grid)
library(tibble)
library(rlang)
library(scales)
library(gridGraphics)

source("~/Desktop/ServerView/MYC/scripts/data_290917/Figures/BindingClasses_expression/ggupset/R/data.R")
source("~/Desktop/ServerView/MYC/scripts/data_290917/Figures/BindingClasses_expression/ggupset/R/scale_upset.R")
source("~/Desktop/ServerView/MYC/scripts/data_290917/Figures/BindingClasses_expression/ggupset/R/theme_combmatrix.R")
source("~/Desktop/ServerView/MYC/scripts/data_290917/Figures/BindingClasses_expression/ggupset/R/axis_combmatrix.R")

# theme_vt2 <- function(){
#   theme_minimal(base_size=18) %+replace%
#     theme(axis.ticks.length = unit(0.3, "cm"),
#           #panel.border = element_rect(size=1, color="black", fill=NA),
#           axis.text = element_text(color="black", face = "bold"),
#           axis.title = element_text(color="black", face = "bold"),
#           legend.text = element_text(color="black", face = "bold"),
#           legend.title = element_text(color="black", face = "bold")
#     )
# }

theme_vt <- function(){
  theme_bw(base_size=18) %+replace%
    theme(axis.ticks.length = unit(0.3, "cm"),
          panel.border = element_rect(size=1.3, color="black", fill=NA),
          axis.text = element_text(color="black"),
          axis.title = element_text(color="black"),
          legend.text = element_text(color="black"),
          legend.title = element_text(color="black")
    )
}

theme_vt2 <- function(){
  theme_minimal(base_size=18) %+replace%
    theme(axis.ticks.length = unit(0.3, "cm"),
          #panel.border = element_rect(size=1, color="black", fill=NA),
          axis.text = element_text(color="black"),
          axis.title = element_text(color="black"),
          legend.text = element_text(color="black"),
          legend.title = element_text(color="black")
    )
}

theme_vt3 <- function(){
  theme_classic(base_size=18) %+replace%
    theme(axis.ticks.length = unit(0.3, "cm"),
          #panel.border = element_rect(size=1, color="black", fill=NA),
          axis.text = element_text(color="black"),
          axis.title = element_text(color="black"),
          legend.text = element_text(color="black"),
          legend.title = element_text(color="black")
    )
}

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

txdb.cstm = readRDS("~/Desktop/ServerView/MYC/results/data_290917/gencode-transcripts/Gencode.v19.ProteinCoding.gr.rds")

OUTPUT = "../figures"

all_bed_lst <- readRDS("../data/Factor_Combinations_beds.rds")

peakAnnoList_txdb <- lapply(all_bed_lst, annotatePeak, TxDb=txdb,
                            tssRegion=c(-1500, 500), verbose=TRUE)

## make a custom plot for ordering in feature plot...
lst_names = names(peakAnnoList_txdb)

annostats_lst = NULL

for(i in 1:length(lst_names)){
  
  my_prot = lst_names[i]
  
  my_df = peakAnnoList_txdb[[my_prot]]@annoStat
  
  annostats_lst[[my_prot]] = my_df
  
}

annostats_lst_df = annostats_lst %>% plyr::ldply() %>%
  dplyr::rename(binding_sites = ".id")

annostats_lst_df2 = annostats_lst_df

## feature_plot2 with upset annotations
annostats_lst_df3 <- annostats_lst_df2

annostats_lst_df3$binding_sites <- factor(annostats_lst_df3$binding_sites, levels = c("MYC_HDAC2_H3K27ac", "MYC_HDAC2", "MYC_H3K27ac",
                                                                                      "HDAC2_H3K27ac", "MYC", "HDAC2", "H3K27ac"))

factor_col <- c("#78CA20", "#AD0AFD", "#F22C1E")
names(factor_col) <- c("MYC", "HDAC2", "H3K27ac")

feature_plot3_upsetAnno <- ggplot(annostats_lst_df3, aes(binding_sites, Frequency, fill = forcats::fct_rev(Feature)))+
  geom_bar(stat="identity", color = "white")+
  scale_fill_brewer(palette="Spectral")+
  theme_vt2()+
  labs(fill='Feature')+
  axis_combmatrix(levels = c("H3K27ac", "HDAC2", "MYC"), expand = TRUE)+
  theme_combmatrix(combmatrix.label.text = element_text(size = 18, color = "black"),
                   combmatrix.label.make_space = TRUE,
                   combmatrix.label.total_extra_spacing = unit(20, "pt"),
                   combmatrix.panel.margin = unit(c(0, 1.5), "pt"),
                   combmatrix.panel.point.size = 6,
                   combmatrix.panel.line.size = 2,
                   combmatrix.panel.point.color.empty = NA)+
  theme(axis.title = element_blank())+
  rotate_y_text(angle = 90, hjust = 1)+
  scale_y_continuous(position = "right")


write.table(annostats_lst_df3, file = "~/Desktop/TRIPLET_FeaturePercentages.txt", quote = F, sep = "\t", row.names = FALSE)

feature_plot3_upsetAnno2 <- ggplot(annostats_lst_df3, aes(binding_sites, Frequency, fill = forcats::fct_rev(Feature)))+
  geom_bar(stat="identity", color = "white")+
  scale_fill_brewer(palette="Spectral")+
  theme_vt3()+
  labs(fill='Feature', y = "% of peaks in feature")+
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), 
                     labels = c("0%", "25%", "50%", "75%", "100%"),
                     expand = expand_scale(mult = c(0, .3)))+
  scale_x_discrete(expand = expand_scale(add = 0.1))+
  axis_combmatrix(levels = c("MYC", "HDAC2", "H3K27ac"), expand = TRUE)+
  theme_combmatrix(combmatrix.label.text = element_text(size = 12, color = "black"),
                   combmatrix.label.make_space = TRUE,
                   combmatrix.label.total_extra_spacing = unit(20, "pt"),
                   combmatrix.panel.margin = unit(c(0, 1.5), "pt"),
                   combmatrix.panel.point.size = 6,
                   combmatrix.panel.line.size = 2,
                   combmatrix.panel.point.color.empty = NA)+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.line.x.bottom = element_blank())

pdf(paste0(OUTPUT, "/MB_Primary_binding_classes_feature_plot_ROTATED_noFrame_upsetAnno.pdf"), height = 10, width = 7.5)
print(feature_plot3_upsetAnno)
dev.off()

save_plot("../figures/MB_Primary_binding_classes_feature_plot_ROTATED_noFrame_upsetAnno.pdf", feature_plot3_upsetAnno2, base_aspect_ratio = 1.3, base_height = 6, base_width = 7)


#plot_grid(upset_recorded, feature_plot3_upsetAnno2, labels = "AUTO", ncol = 1, align = "h")

# upset_recorded + feature_plot3_upsetAnno2 + patchwork::plot_layout(ncol = 2)
