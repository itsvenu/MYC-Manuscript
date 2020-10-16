## source: /b06x-lsdf/infra5-lsdf/mbHDAC2/scripts/data_290917/Figures/BindingClasses_expression/bindingClasses_expression_analysis.R
## different binding classes - gep

setwd("~/Desktop/ServerView/MYC/scripts/MycHdaci_code/scripts")

library(tidyverse)
library(data.table)
library(ggpubr)
library(EnvStats)

library(gtable)
library(grid)
library(tibble)
library(rlang)
library(scales)

source("~/Desktop/ServerView/MYC/scripts/data_290917/Figures/BindingClasses_expression/ggupset/R/data.R")
source("~/Desktop/ServerView/MYC/scripts/data_290917/Figures/BindingClasses_expression/ggupset/R/scale_upset.R")
source("~/Desktop/ServerView/MYC/scripts/data_290917/Figures/BindingClasses_expression/ggupset/R/theme_combmatrix.R")
source("~/Desktop/ServerView/MYC/scripts/data_290917/Figures/BindingClasses_expression/ggupset/R/axis_combmatrix.R")

# theme_vt <- function(){
#   theme_bw(base_size=18) %+replace%
#     theme(axis.ticks.length = unit(0.3, "cm"),
#           panel.border = element_rect(size=1, color="black", fill=NA),
#           axis.text = element_text(color="black", face = "bold"),
#           axis.title = element_text(color="black", face = "bold"),
#           legend.text = element_text(color="black", face = "bold"),
#           legend.title = element_text(color="black", face = "bold")
#     )
# }

theme_vt <- function(){
  theme_bw(base_size=18) %+replace%
    theme(axis.ticks.length = unit(0.3, "cm"),
          panel.border = element_rect(size=1, color="black", fill=NA),
          axis.text = element_text(color="black"),
          axis.title = element_text(color="black"),
          legend.text = element_text(color="black"),
          legend.title = element_text(color="black")
    )
}

INPUT = "~/Desktop/ServerView/MYC/results/data_290917/Figures/BindingClasses_expression/gene_sets"

OUTPUT = "../figures"


matched_exp_tpm <- readRDS("../data/Primary_ChIP_matched_RNAseq_TPM.rds")

matched_pid_mean <- matched_exp_tpm %>% 
  reshape2::melt() %>% 
  dplyr::mutate(log_tpm = log2(value+1)) %>% 
  dplyr::group_by(target_id) %>% 
  dplyr::summarise(MEAN = mean(log_tpm)) %>% 
  as.data.frame()

## each file contains list of genes symbols on per line
## file name contains which binding classes each gene set belong to.

all_gene_sets = readRDS("../data/Primary_BindingClasses_genesets_paths.rds")

all_gene_sets_matched_mlt = NULL

for(i in 1:length(all_gene_sets)){
  
  d = all_gene_sets[i]
  
  d_name = gsub("_promoter_overlaping_uniqgenes.txt", "", basename(d))
  
  d_dat = read.delim(d, header = FALSE)
  
  d_class_exp = merge(matched_pid_mean, d_dat, by.x = "target_id", by.y = "V1") %>%
    dplyr::mutate(CLASS = d_name)
  
  all_gene_sets_matched_mlt = rbind(all_gene_sets_matched_mlt, d_class_exp)
  
}  

## we need NONE gene set
## read all protein coding genes
## extract those which are not overlapping with any class
gencode_genes = read.delim("~/Desktop/ServerView/MYC/results/data_290917/gencode-transcripts/Gencode.v19.ProteinCoding.genelist.txt", header = FALSE)

all_gene_classes_exp = setdiff(gencode_genes$V1, all_gene_sets_matched_mlt$target_id) %>%
  as.data.frame() %>%
  dplyr::rename(NONE = ".") %>%
  merge(matched_pid_mean, by.x = "NONE", by.y = "target_id") %>%
  dplyr::mutate(CLASS = "NONE") %>%
  dplyr::rename(target_id = NONE) %>%
  rbind(all_gene_sets_matched_mlt)

saveRDS(all_gene_classes_exp, file = "../data/Primary_AllClasses_matched_pid_exp.rds")

## fix x-axis order..

all_gene_classes_exp$CLASS = factor(all_gene_classes_exp$CLASS,
                                    levels = c("MYC_HDAC2_H3K27ac", "MYC_HDAC2", "MYC_H3K27ac",
                                               "HDAC2_H3K27ac", "MYC", "HDAC2", "H3K27ac", "NONE"))

factor_col <- c("#78CA20", "#AD0AFD", "#F22C1E", "lightgray")
names(factor_col) <- c("MYC", "HDAC2", "H3K27ac", "NONE")

bindingClasses_exp_plt <- all_gene_classes_exp %>%
  dplyr::filter(MEAN > 0.5) %>%
  ggboxplot(x = "CLASS", y = "MEAN", fill = "CLASS")+
  scale_fill_brewer(palette = "Set1")+
  theme_vt()+
  ylab("Gene expression (log2 TPM)")+
  theme(axis.title.x = element_blank(),
        legend.position = "bottom")+
  stat_compare_means(ref.group = "MYC_HDAC2_H3K27ac", label = "p.signif", size = 7, color = "red")+
  rotate_x_text(angle = 35, hjust = 1)+
  stat_n_text(size = 6, fontface = "bold", text.box = TRUE)+
  axis_combmatrix(levels = c("MYC", "HDAC2", "H3K27ac", "NONE"))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())+
  theme_combmatrix(combmatrix.label.text = element_text(size = 18, color = "black"),
                   combmatrix.label.make_space = FALSE,
                   combmatrix.label.total_extra_spacing = unit(20, "pt"),
                   combmatrix.panel.margin = unit(c(0.005, 1.5), "pt"),
                   combmatrix.panel.point.size = 6,
                   combmatrix.panel.line.size = 2,
                   combmatrix.panel.point.color.empty = NA)

## change theme
bindingClasses_exp_plt2 <- all_gene_classes_exp %>%
  dplyr::filter(MEAN > 0.5) %>%
  ggboxplot(x = "CLASS", y = "MEAN", fill = "CLASS")+
  scale_fill_brewer(palette = "Set1")+
  theme_classic(base_size = 18)+
  ylab("Gene expression (Log2 TPM)")+
  theme(axis.title.x = element_blank(),
        legend.position = "bottom")+
  stat_compare_means(ref.group = "MYC_HDAC2_H3K27ac", label = "p.signif", size = 7, color = "red")+
  rotate_x_text(angle = 35, hjust = 1)+
  stat_n_text(size = 4, text.box = TRUE)+
  scale_x_discrete(expand = expand_scale(add = 0.3))+
  axis_combmatrix(levels = c("MYC", "HDAC2", "H3K27ac", "NONE"))+
  theme(axis.text.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.major.y = element_line(size = 0.8, color = "lightgray"),
        axis.ticks.length = unit(0.3, "cm"))+
  theme_combmatrix(combmatrix.label.text = element_text(size = 18, color = "black"),
                   combmatrix.label.make_space = FALSE,
                   combmatrix.label.total_extra_spacing = unit(20, "pt"),
                   combmatrix.panel.margin = unit(c(0.005, 1.5), "pt"),
                   combmatrix.panel.point.size = 6,
                   combmatrix.panel.line.size = 2,
                   combmatrix.panel.point.color.empty = NA)

## proportional bar plot of classes
all_genes_donut <- all_gene_classes_exp %>% dplyr::count(CLASS) %>%
  dplyr::rename(all_genes = n)

## add expressed genes
all_genes_donut <- all_gene_classes_exp %>%
  dplyr::filter(MEAN > 0.5) %>%
  dplyr::count(CLASS) %>%
  dplyr::rename(expressed_genes = n) %>%
  merge(all_genes_donut, by = "CLASS") %>%
  dplyr::mutate(nonexpressed_genes = all_genes - expressed_genes)

## use this data to plot donut plot using python
all_genes_donut_mlt <- all_genes_donut %>%
  dplyr::select(-c(all_genes)) %>%
  reshape2::melt() %>%
  dplyr::rename(parent = CLASS, node = variable, size = value)

##
exp_prop <- all_genes_donut_mlt %>%
  dplyr::filter(node == "expressed_genes") %>%
  dplyr::mutate(prop = size/sum(size))

nonexp_prop <- all_genes_donut_mlt %>%
  dplyr::filter(node == "nonexpressed_genes") %>%
  dplyr::mutate(prop = size/sum(size))

# > all_genes_donut$expressed_genes %>% sum
# [1] 6444
# > all_genes_donut$nonexpressed_genes %>% sum
# [1] 13785

all_prop <- rbind(exp_prop, nonexp_prop)
all_prop <- all_prop %>%
  dplyr::mutate(exp_class = case_when(node == "expressed_genes" ~ "Active\n(n=6,444)",
                                      node == "nonexpressed_genes" ~ "Silent\n(n=13,785)"))

all_prop$exp_class <- factor(all_prop$exp_class, levels = c("Silent\n(n=13,785)", "Active\n(n=6,444)"))
all_prop$parent <- factor(all_prop$parent, levels = c("MYC_HDAC2_H3K27ac", "MYC_HDAC2",
                                                      "MYC_H3K27ac", "HDAC2_H3K27ac", "MYC", "HDAC2", "H3K27ac", "NONE"))

sb_cols <- RColorBrewer::brewer.pal(n = 8, "Set1")
names(sb_cols) <- all_genes_donut$CLASS %>% levels()

active_silent_propBar <- ggbarplot(all_prop, x = "exp_class", y = "prop", fill = "parent", palette = sb_cols, color = NA,
          alpha = "exp_class")+
  scale_alpha_discrete(range = c(0.4, 1))+
  xlab("")+ylab("% of protein coding genes")+
  theme_classic(base_size = 18)+
  theme(axis.text = element_text(color="black"),
        axis.title = element_text(color="black"),
        legend.text = element_text(color="black"),
        legend.title = element_blank(),
        axis.ticks.length = unit(0.3, "cm"))+
  guides(alpha = FALSE)+
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), 
                     labels = c("0%", "25%", "50%", "75%", "100%"))

active_silent_propBar2 <- ggbarplot(all_prop, x = "exp_class", y = "prop", 
                                    fill = "parent", palette = sb_cols, color = NA,
                                    alpha = "exp_class", width = 0.9)+
  scale_alpha_discrete(range = c(0.4, 1))+
  xlab("")+ylab("% of protein coding genes")+
  theme_classic(base_size = 18)+
  theme(axis.text = element_text(color="black", size = 12),
        axis.title = element_text(color="black", size = 14),
        legend.text = element_text(color="black"),
        legend.title = element_blank(),
        axis.line.x.bottom = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        axis.ticks.x = element_blank())+
  guides(alpha = FALSE)+
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1),
                     labels = c("0%", "25%", "50%", "75%", "100%"),
                     expand = expand_scale(mult = c(0, 0)))+
  scale_x_discrete(expand = expand_scale(add = 0.5))

## both in one
pdf("../figures/Primary_BindingClasses_exp_BoxPropBar.pdf", height = 7.5, width = 12)

active_silent_propBar +
  theme(legend.position = "bottom")+
  guides(fill=guide_legend(ncol=2,byrow=FALSE)) +
  bindingClasses_exp_plt + 
  patchwork::plot_layout(widths = c(1, 3))

dev.off()

## single plots
pdf("../figures/Primary_BindingClasses_Active_vs_Silent_proportionalBarPlot.pdf", height = 7, width = 4.5)

active_silent_propBar +
  theme(legend.position = "none")

active_silent_propBar +
  theme(legend.position = "bottom")+
  guides(fill=guide_legend(ncol=2,byrow=FALSE))

dev.off()

## boxplots
pdf("../figures/Primary_BindingClasses_ActiveGenes_exp_boxplot.pdf", height = 7, width = 11)

bindingClasses_exp_plt

dev.off()

###
save_plot("../figures/Primary_BindingClasses_Active_vs_Silent_proportionalBarPlot_1.pdf", active_silent_propBar2, base_aspect_ratio = 1.3, base_height = 5)
save_plot("../figures/Primary_BindingClasses_Active_vs_Silent_proportionalBarPlot_2.pdf", active_silent_propBar2+theme(legend.position = "none"), base_aspect_ratio = 1.3, base_height = 5, base_width = 4)

save_plot("../figures/Primary_BindingClasses_ActiveGenes_exp_boxplot.pdf", bindingClasses_exp_plt, base_aspect_ratio = 1.3, base_height = 6)
save_plot("../figures/Primary_BindingClasses_ActiveGenes_exp_boxplot_classic.pdf", bindingClasses_exp_plt2, base_aspect_ratio = 1.3, base_height = 6)

## only triplet bound - expressed vs non-expressed
triplet_all <- all_gene_classes_exp %>% 
  dplyr::mutate(exp = dplyr::case_when(MEAN <= 0.5 ~ "Silent",
                                       MEAN > 0.5 ~ "Active")) %>% 
  dplyr::filter(CLASS == "MYC_HDAC2_H3K27ac")

# > dim(triplet_all)
# [1] 10117     4
triplet_all <- triplet_all %>% 
  dplyr::count(exp) %>% 
  dplyr::mutate(fraction = n/10117) %>% 
  dplyr::mutate(group = "Triplet all")

triplet_all$lab <- c("Active (47.1%)", "Silenet (52.9%)")

# A tibble: 2 x 3
# exp        n fraction
# <chr>  <int>    <dbl>
# 1 Active  4767    0.471
# 2 Silent  5350    0.529

pdf("../figures/Primary_TripletAll_ActiveSilet_pie.pdf", height = 5, width = 6)

ggpubr::ggpie(triplet_all, x = "fraction", 
              label = "lab", 
              lab.pos = "in",
              lab.font = c(7, "bold", "black"),
              fill = "exp", 
              color = "white",
              palette = c("Active" = "indianred", "Silent" = "gray"))+
  theme(legend.position = "none")

dev.off()

###
## triplet - doublet - singlet - none comparison
library(gginnards)

trip_dub_sings <- all_gene_classes_exp %>% 
  dplyr::mutate(group = dplyr::case_when(CLASS == "MYC_HDAC2_H3K27ac" ~ "Triplet",
                                         CLASS == "MYC_HDAC2" | CLASS == "MYC_H3K27ac" | CLASS == "HDAC2_H3K27ac" ~ "Doublet",
                                         CLASS == "MYC" | CLASS == "HDAC2" | CLASS == "H3K27ac" ~ "Singlet",
                                         CLASS == "NONE" ~ "NONE")) %>% 
  dplyr::filter(MEAN > 0.5) 

trip_dub_sings$group <- factor(trip_dub_sings$group, levels = c("Triplet", "Doublet", "Singlet", "NONE"))

my_comp <- list(c("Triplet", "Doublet"),
                c("Triplet", "Singlet"),
                c("Triplet", "NONE"),
                c("Doublet", "Singlet"),
                c("Doublet", "NONE"))

trip_dub_sings_plt <- trip_dub_sings %>% 
  ggboxplot(x = "group", y = "MEAN", fill = "group")+
  stat_compare_means(comparisons = my_comp, size = 0.9)+
  scale_fill_manual(values = c("Triplet" = "#A93226", "Doublet" = "#CD6155", "Singlet" = "#E6B0AA", "NONE" = "gray"))+
  theme_classic(base_size = 18)+
  theme(axis.ticks.length = unit(0.3, "cm"),
        axis.text = element_text(color="black"),
        legend.position = "none")+
  EnvStats::stat_n_text(text.box = TRUE, size = 6)+
  ylab("Gene expression (Log2 TPM)")+
  xlab("Binding status of genes")

trip_dub_sings_plt$layers[[2]]$aes_params$textsize <- 6

pdf("../figures/Primary_Triplet_Doublet_Singlet_NONE_expComparison.pdf", height = 6, width = 7)
print(trip_dub_sings_plt)
dev.off()


