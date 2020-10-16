## source: /b06x-lsdf/infra5-lsdf/mbHDAC2/scripts/affy-ge/DMSO_HIGHDOSE_VOLCANO_theme_vt.R
## source: /b06x-lsdf/infra5-lsdf/mbHDAC2/scripts/affy-ge/all_timepoints_deg/closestEnhancers/metaplot_smoothing.R
## source: /b06x-lsdf/infra5-lsdf/mbHDAC2/results/affy-ge/all_timepoints_deg/all_timepoints_deg_coords/AllGenes_TES_RNApolII/allGenes_TES_RNApolII.R

## volcano plot of GEP with highlighting binding status from ChIP
## status of DEG interms of protein binding before treatment

setwd("~/Desktop/ServerView/MYC/scripts/MycHdaci_code/scripts")

library(tidyverse)
library(EnvStats)
library(ggpubr)

theme_vt <- function(){
  theme_bw(base_size=18) %+replace%
    theme(axis.ticks.length = unit(0.3, "cm"),
          panel.border = element_rect(size=1.3, color="black", fill=NA),
          axis.text = element_text(color="black", face = "bold"),
          axis.title = element_text(color="black", face = "bold"),
          legend.text = element_text(color="black", face = "bold"),
          legend.title = element_text(color="black", face = "bold")
    )
}
theme_vt_ms <- function(){
  theme_bw(base_size=18) %+replace%
    theme(axis.ticks.length = unit(0.3, "cm"),
          panel.border = element_rect(size=1.75, color="black", fill=NA),
          axis.text = element_text(color="black"),
          axis.title = element_text(color="black"),
          legend.text = element_text(color="black"),
          legend.title = element_text(color="black")
    )
}

dat = readRDS("../data/DMSO_MS275_6h_GEP_limmaResults.rds")

dat_fmt2 = dat %>% dplyr::mutate(sig=ifelse(dat$adj.P.Val<0.1 & abs(dat$logFC) > 1, "FDR<0.1", "FDR>0.1")) 

## triplet genes from HD-MB03
cline_triplet = read.delim("~/Desktop/ServerView/MYC/results/data_30052018/downstreamAnalysis/treated_untreated_intervene/Untreated_vs_treated/analysis_12092018_ds/downstreamAnalysis/analysis_12102018/04_geneAssignment/assignedGenes/unique_sets/uniq_sets_newsets/UNIQ_SETS_NEWSETS_07122018/classes/1kb_quant/RANKING_SETS/MYC_HDAC2_H3K27ac_ALL_GENES_wholegeneCoordinates.txt", header = FALSE) %>% 
  dplyr::select(V4) %>% unique()

cline_triplet_dat_fmt2 = intersect(cline_triplet$V4, dat_fmt2$gene.symbols) %>% 
  as.data.frame() %>% 
  dplyr::rename(GENE = ".") %>% 
  merge(dat_fmt2, by.x = "GENE", by.y = "gene.symbols") %>% 
  dplyr::mutate(TRIPLET = " & bound by triplet")

cline_NOTtriplet_dat_fmt2 = setdiff(dat_fmt2$gene.symbols, cline_triplet$V4) %>% 
  as.data.frame() %>% 
  dplyr::rename(GENE = ".") %>% 
  merge(dat_fmt2, by.x = "GENE", by.y = "gene.symbols") %>% 
  dplyr::mutate(TRIPLET = " & not bound by triplet")

## merge both into one df
dat_fmt2_clineTriplet = rbind(cline_triplet_dat_fmt2, cline_NOTtriplet_dat_fmt2) %>% 
  dplyr::mutate(sig_triplet = paste0(sig, TRIPLET))

cline_vlcno_clr_alltriplet_colord_alphd = ggplot(dat_fmt2_clineTriplet, aes(logFC, -log10(adj.P.Val), color = sig_triplet, alpha = sig_triplet)) +
  geom_point() +
  scale_color_manual(values=c("FDR>0.1 & bound by triplet" = "indianred", "FDR<0.1 & bound by triplet" = "red", "FDR>0.1 & not bound by triplet" = "gray", "FDR<0.1 & not bound by triplet" = "black"))+
  scale_alpha_discrete(range=rev(c(0.4, 1)))+
  geom_hline(yintercept = 1, linetype = 2)+
  geom_vline(xintercept = c(1, -1), linetype = 2)+
  theme_vt()+
  xlab("Log2 fold change")+
  ylab("-log10(FDR)")+
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom")+
  ggtitle("up-442, down-316; 6h; HD-MB03")+
  guides(color=guide_legend(nrow=4,byrow=TRUE))

cline_vlcno_clr_alltriplet_colord_alphd_ms = ggplot(dat_fmt2_clineTriplet, aes(logFC, -log10(adj.P.Val), color = sig_triplet, alpha = sig_triplet)) +
  geom_point() +
  scale_color_manual(values=c("FDR>0.1 & bound by triplet" = "indianred", 
                              "FDR<0.1 & bound by triplet" = "red", 
                              "FDR>0.1 & not bound by triplet" = "gray", 
                              "FDR<0.1 & not bound by triplet" = "black"))+
  scale_alpha_discrete(range=rev(c(0.4, 1)))+
  geom_hline(yintercept = 1, linetype = 2)+
  geom_vline(xintercept = c(1, -1), linetype = 2)+
  theme_vt_ms()+
  xlab("Log2 fold change")+
  ylab("-log10(FDR)")+
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom")+
  ggtitle("up-442, down-316; 6h; HD-MB03")+
  guides(color=guide_legend(nrow=4,byrow=TRUE))

## change colors

pdf("../figures/DMSO_MS275_6h_GEP_limmaResults_volcanoPlot_tripletHighlighted_colorChange.pdf", height = 7, width = 5)

ggplot(dat_fmt2_clineTriplet, aes(logFC, -log10(adj.P.Val), color = sig_triplet, alpha = sig_triplet)) +
  geom_point(size = 0.8) +
  scale_color_manual(values=c("FDR>0.1 & bound by triplet" = "indianred", 
                              "FDR<0.1 & bound by triplet" = "red", 
                              "FDR>0.1 & not bound by triplet" = "gray", 
                              "FDR<0.1 & not bound by triplet" = "#808000"))+
  scale_alpha_discrete(range=rev(c(0.4, 1)))+
  geom_hline(yintercept = 1, linetype = 2)+
  geom_vline(xintercept = c(1, -1), linetype = 2)+
  theme_vt_ms()+
  xlab("Log2 fold change")+
  ylab("-log10(FDR)")+
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom")+
  ggtitle("up-442, down-316; 6h; HD-MB03")+
  guides(color=guide_legend(nrow=4,byrow=TRUE))

dev.off()

### update volcano plot
dat_fmt2_clineTriplet2 <- dat_fmt2_clineTriplet %>% 
  dplyr::mutate(sig_triplet_fc = case_when(logFC < 0 ~ paste0(sig_triplet, " & negtive_sig"),
                                           logFC > 0 ~ paste0(sig_triplet, " & positive_sig")))

pdf("../figures/DMSO_MS275_6h_GEP_limmaResults_volcanoPlot_tripletHighlighted_colorChangesBluesReds.pdf", height = 10, width = 7)

ggplot(dat_fmt2_clineTriplet2, aes(logFC, -log10(adj.P.Val), 
                                   color = sig_triplet_fc, 
                                   alpha = sig_triplet_fc)) +
  geom_point(size = 0.8) +
  scale_color_manual(values=c("FDR>0.1 & bound by triplet & positive_sig" = "indianred",
                              "FDR>0.1 & bound by triplet & negtive_sig" = "lightblue",
                              "FDR<0.1 & bound by triplet & negtive_sig" = "blue",
                              "FDR<0.1 & bound by triplet & positive_sig" = "red",
                              "FDR>0.1 & not bound by triplet & positive_sig" = "black",
                              "FDR>0.1 & not bound by triplet & negtive_sig" = "black",
                              "FDR<0.1 & not bound by triplet & positive_sig" = "black",
                              "FDR<0.1 & not bound by triplet & negtive_sig" = "black"))+
  scale_alpha_discrete(range=rev(c(0.4, 1)))+
  geom_hline(yintercept = 1, linetype = 2)+
  geom_vline(xintercept = c(1, -1), linetype = 2)+
  theme_vt_ms()+
  xlab("Log2 fold change")+
  ylab("-log10(FDR)")+
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom")+
  #ggtitle("up-442, down-316; 6h; HD-MB03")+
  guides(color=guide_legend(ncol = 1))

ggplot(dat_fmt2_clineTriplet2, aes(logFC, -log10(adj.P.Val), 
                                   color = sig_triplet_fc, 
                                   alpha = sig_triplet_fc)) +
  geom_point(size = 0.8) +
  scale_color_manual(values=c("FDR>0.1 & bound by triplet & positive_sig" = "indianred",
                              "FDR>0.1 & bound by triplet & negtive_sig" = "lightblue",
                              "FDR<0.1 & bound by triplet & negtive_sig" = "blue",
                              "FDR<0.1 & bound by triplet & positive_sig" = "red",
                              "FDR>0.1 & not bound by triplet & positive_sig" = "darkgray",
                              "FDR>0.1 & not bound by triplet & negtive_sig" = "darkgray",
                              "FDR<0.1 & not bound by triplet & positive_sig" = "darkgray",
                              "FDR<0.1 & not bound by triplet & negtive_sig" = "darkgray"))+
  scale_alpha_discrete(range=rev(c(0.4, 1)))+
  geom_hline(yintercept = 1, linetype = 2)+
  geom_vline(xintercept = c(1, -1), linetype = 2)+
  theme_vt_ms()+
  xlab("Log2 fold change")+
  ylab("-log10(FDR)")+
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom")+
  #ggtitle("up-442, down-316; 6h; HD-MB03")+
  guides(color=guide_legend(ncol = 1))

ggplot(dat_fmt2_clineTriplet2, aes(logFC, -log10(adj.P.Val), 
                                   color = sig_triplet_fc, 
                                   alpha = sig_triplet_fc)) +
  geom_point(size = 0.8) +
  scale_color_manual(values=c("FDR>0.1 & bound by triplet & positive_sig" = "indianred",
                              "FDR>0.1 & bound by triplet & negtive_sig" = "lightblue",
                              "FDR<0.1 & bound by triplet & negtive_sig" = "blue",
                              "FDR<0.1 & bound by triplet & positive_sig" = "red",
                              "FDR>0.1 & not bound by triplet & positive_sig" = "darkgray",
                              "FDR>0.1 & not bound by triplet & negtive_sig" = "darkgray",
                              "FDR<0.1 & not bound by triplet & positive_sig" = "black",
                              "FDR<0.1 & not bound by triplet & negtive_sig" = "black"))+
  scale_alpha_discrete(range=rev(c(0.4, 1)))+
  geom_hline(yintercept = 1, linetype = 2)+
  geom_vline(xintercept = c(1, -1), linetype = 2)+
  theme_vt_ms()+
  xlab("Log2 fold change")+
  ylab("-log10(FDR)")+
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom")+
  #ggtitle("up-442, down-316; 6h; HD-MB03")+
  guides(color=guide_legend(ncol = 1))

ggplot(dat_fmt2_clineTriplet2, aes(logFC, -log10(adj.P.Val), 
                                   color = sig_triplet_fc)) +
  geom_point(size = 0.8) +
  scale_color_manual(values=c("FDR>0.1 & bound by triplet & positive_sig" = "indianred",
                              "FDR>0.1 & bound by triplet & negtive_sig" = "lightblue",
                              "FDR<0.1 & bound by triplet & negtive_sig" = "blue",
                              "FDR<0.1 & bound by triplet & positive_sig" = "red",
                              "FDR>0.1 & not bound by triplet & positive_sig" = "gray80",
                              "FDR>0.1 & not bound by triplet & negtive_sig" = "gray80",
                              "FDR<0.1 & not bound by triplet & positive_sig" = "black",
                              "FDR<0.1 & not bound by triplet & negtive_sig" = "black"))+
  #scale_alpha_discrete(range=rev(c(0.4, 1)))+
  geom_hline(yintercept = 1, linetype = 2)+
  geom_vline(xintercept = c(1, -1), linetype = 2)+
  theme_vt_ms()+
  xlab("Log2 fold change")+
  ylab("-log10(FDR)")+
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom")+
  #ggtitle("up-442, down-316; 6h; HD-MB03")+
  guides(color=guide_legend(ncol = 1))

dev.off()


# ggplot(dat_fmt2_clineTriplet, aes(logFC, -log10(adj.P.Val), color = sig_triplet, alpha = sig_triplet)) +
#   geom_point() +
#   scale_color_manual(values=c("FDR>0.1 & bound by triplet" = "indianred", 
#                               "FDR<0.1 & bound by triplet" = "red", 
#                               "FDR>0.1 & not bound by triplet" = "gray", 
#                               "FDR<0.1 & not bound by triplet" = "lightgreen"))+
#   scale_alpha_discrete(range=rev(c(0.4, 1)))+
#   geom_hline(yintercept = 1, linetype = 2)+
#   geom_vline(xintercept = c(1, -1), linetype = 2)+
#   theme_vt_ms()+
#   xlab("Log2 fold change")+
#   ylab("-log10(FDR)")+
#   theme(legend.title = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         legend.position = "bottom")+
#   ggtitle("up-442, down-316; 6h; HD-MB03")+
#   guides(color=guide_legend(nrow=4,byrow=TRUE))

pdf("../figures/DMSO_MS275_6h_GEP_limmaResults_volcanoPlot_tripletHighlighted.pdf", height = 9, width = 6.5)
print(cline_vlcno_clr_alltriplet_colord_alphd)
print(cline_vlcno_clr_alltriplet_colord_alphd_ms)
dev.off()

## DEG protein binding status before treatment
source("~/Desktop/ServerView/MYC/scripts/data_30052018/downstreamAnalysis/Untreated_treated_analysis/analysis_12102018/04_geneAssignment/gene_downstream/uniq_sets_plots/HDACi_lib_modified_functions.R")
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

## overlapping flat violin function
return_faltviolin_dmsoMs275 <- function(dat_df){
  
  # dat_df <- cl_triplet_up %>% 
  #   dplyr::select(GENE) %>% 
  #   merge(all_quant_pergene2, by.x = "GENE", by.y = "V4") %>% 
  #   dplyr::select(-c(matches("RNA")))
  dat_df_mlt <- dat_df %>% 
    reshape2::melt() %>% 
    dplyr::mutate(protein = variable) %>% 
    dplyr::mutate(protein = gsub("_.*", "", protein)) 
  
  ## fix some orders
  dat_df_mlt$protein <- factor(dat_df_mlt$protein, levels = c("MYC", "HDAC2", "H3K27ac"))
  dat_df_mlt$variable <- factor(dat_df_mlt$variable, levels = c("MYC_untreated", "MYC_treated",
                                                                "HDAC2_untreated", "HDAC2_treated",
                                                                "H3K27ac_untreated", "H3K27ac_treated"))
  
  plt <- ggplot(data = dat_df_mlt, aes(y = value, x = protein, fill = variable)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .6)+
    #geom_point(aes(y = value, color = DEG), position = position_jitter(width = .05), size = 1, alpha = 0.8) +
    geom_boxplot(width = .3, outlier.shape = NA) +
    guides(fill = FALSE) +
    guides(color = FALSE)+
    scale_fill_manual(values = c("MYC_untreated" = "darkgray", "MYC_treated" = "#78CA20",
                                 "HDAC2_untreated" = "darkgray", "HDAC2_treated" = "#AD0AFD",
                                 "H3K27ac_untreated" = "darkgray", "H3K27ac_treated" = "#F22C1E"))+
    scale_color_manual(values = c("MYC_untreated" = "darkgray", "MYC_treated" = "#78CA20",
                                  "HDAC2_untreated" = "darkgray", "HDAC2_treated" = "#AD0AFD",
                                  "H3K27ac_untreated" = "darkgray", "H3K27ac_treated" = "#F22C1E"))+
    theme_vt()+
    ylab("input normalized IP signal (log2)")+
    theme(axis.title.x = element_blank(),
          legend.title = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())+
    stat_compare_means(aes(group = variable), label = "p.format", size = 6, label.y = max(dat_df_mlt$value) + 0.5)
  
  return(plt)
  
  
}

## up
cl_all_up = dat_fmt2_clineTriplet %>% dplyr::filter(sig != "FDR>0.1") %>% dplyr::filter(logFC > 1) %>% 
  dplyr::select(GENE, TRIPLET) %>% unique()

cl_triplet_up = dat_fmt2_clineTriplet %>% dplyr::filter(sig != "FDR>0.1") %>% dplyr::filter(logFC > 1) %>% 
  dplyr::select(GENE, TRIPLET) %>% unique() %>% 
  dplyr::filter(TRIPLET == " & bound by triplet")

## down
cl_all_down = dat_fmt2_clineTriplet %>% dplyr::filter(sig != "FDR>0.1") %>% dplyr::filter(logFC < -1) %>% 
  dplyr::select(GENE, TRIPLET) %>% unique()

cl_triplet_down = dat_fmt2_clineTriplet %>% dplyr::filter(sig != "FDR>0.1") %>% dplyr::filter(logFC < -1) %>% 
  dplyr::select(GENE, TRIPLET) %>% unique() %>% 
  dplyr::filter(TRIPLET == " & bound by triplet")

all_quant_pergene2 = readRDS("../data/HDMB03_DMSO_MS275_ProteinEnrichments_ProteinCodingPromoters.rds")

# cl_triplet_up %>% 
#   dplyr::select(GENE) %>% 
#   merge(all_quant_pergene2, by.x = "GENE", by.y = "V4") %>% 
#   dplyr::select(-c(matches("RNA"))) %>% 
#   return_faltviolin_dmsoMs275()+
#   ggtitle("MS275-UP, n=299/397")+
#   theme_classic(base_size = 20)+
#   xlab("")+
#   scale_x_discrete(expand = expand_scale(add = 0.35))+
#   theme(axis.text = element_text(color = "black"),
#         axis.ticks.length = unit(0.3, "cm"))


pdf("../figures/HDMB03_6h_GEP_SIG_DEG_BINDINGPATTERNS_onlyTripletBound.pdf", height = 5, width = 8)

cl_triplet_up %>% 
  dplyr::select(GENE) %>% 
  merge(all_quant_pergene2, by.x = "GENE", by.y = "V4") %>% 
  dplyr::select(-c(matches("RNA"))) %>% 
  return_faltviolin_dmsoMs275()+
  ggtitle("MS275-UP, n=299/397")

cl_triplet_down %>% 
  dplyr::select(GENE) %>% 
  merge(all_quant_pergene2, by.x = "GENE", by.y = "V4") %>% 
  dplyr::select(-c(matches("RNA"))) %>%
  return_faltviolin_dmsoMs275()+
  ggtitle("MS275-DOWN, n=262/301")

dev.off()

## classic theme
pdf("../figures/HDMB03_6h_GEP_SIG_DEG_BINDINGPATTERNS_onlyTripletBound_classic.pdf", height = 4.5, width = 6)

cl_triplet_up %>% 
  dplyr::select(GENE) %>% 
  merge(all_quant_pergene2, by.x = "GENE", by.y = "V4") %>% 
  dplyr::select(-c(matches("RNA"))) %>% 
  return_faltviolin_dmsoMs275()+
  ggtitle("MS275-UP, n=299/397")+
  theme_classic(base_size = 20)+
  xlab("")+
  scale_x_discrete(expand = expand_scale(add = 0.35))+
  theme(axis.text = element_text(color = "black"),
        axis.ticks.length = unit(0.3, "cm"))

cl_triplet_down %>% 
  dplyr::select(GENE) %>% 
  merge(all_quant_pergene2, by.x = "GENE", by.y = "V4") %>% 
  dplyr::select(-c(matches("RNA"))) %>%
  return_faltviolin_dmsoMs275()+
  ggtitle("MS275-DOWN, n=262/301")+
  theme_classic(base_size = 20)+
  xlab("")+
  scale_x_discrete(expand = expand_scale(add = 0.35))+
  theme(axis.text = element_text(color = "black"),
        axis.ticks.length = unit(0.3, "cm"))

dev.off()

## plot the differences for 3 factors

pdf("../figures/HDMB03_6h_GEP_SIG_DEG_BINDINGPATTERNS_onlyTripletBound_BindingDifference.pdf", height = 5, width = 4.5)

cl_triplet_up %>% 
  dplyr::select(GENE) %>% 
  merge(all_quant_pergene2, by.x = "GENE", by.y = "V4") %>% 
  dplyr::select(-c(matches("RNA"))) %>% 
  dplyr::mutate(MYC = MYC_treated - MYC_untreated,
                HDAC2 = HDAC2_treated - HDAC2_untreated,
                H3K27ac = H3K27ac_treated - H3K27ac_untreated) %>% 
  dplyr::select(-c(matches("treated"))) %>% 
  reshape2::melt() %>% 
  ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable")+
  scale_fill_manual(values = c("MYC" = "#78CA20", "HDAC2" = "#AD0AFD", "H3K27ac" = "#F22C1E"))+
  theme_classic(base_size = 18)+
  theme(legend.position = "none",
        axis.text = element_text(color="black"),
        axis.ticks.length = unit(0.3, "cm"),
        axis.ticks = element_line(color = "black"))+
  ylab("Binding difference (entinostat-DMSO)")+xlab("")+
  ggtitle("Triplet-up")
  
  
cl_triplet_down %>% 
  dplyr::select(GENE) %>% 
  merge(all_quant_pergene2, by.x = "GENE", by.y = "V4") %>% 
  dplyr::select(-c(matches("RNA"))) %>% 
  dplyr::mutate(MYC = MYC_treated - MYC_untreated,
                HDAC2 = HDAC2_treated - HDAC2_untreated,
                H3K27ac = H3K27ac_treated - H3K27ac_untreated) %>% 
  dplyr::select(-c(matches("treated"))) %>% 
  reshape2::melt() %>% 
  ggpubr::ggboxplot(x = "variable", y = "value", fill = "variable")+
  scale_fill_manual(values = c("MYC" = "#78CA20", "HDAC2" = "#AD0AFD", "H3K27ac" = "#F22C1E"))+
  theme_classic(base_size = 18)+
  theme(legend.position = "none",
        axis.text = element_text(color="black"),
        axis.ticks.length = unit(0.3, "cm"),
        axis.ticks = element_line(color = "black"))+
  ylab("Binding difference (entinostat-DMSO)")+xlab("")+
  ggtitle("Triplet-down")

dev.off()

## non-triplet up/down plots
nonTriplet_up <- read.delim("~/Desktop/ServerView/MYC/results/affy-ge/DEG_DMEA/E_BOX/NONTRIPLET_DEG/NON_TRIPLET_DEG_up_tss.bed", header = FALSE) %>% 
  dplyr::select(V4)

nonTriplet_down <- read.delim("~/Desktop/ServerView/MYC/results/affy-ge/DEG_DMEA/E_BOX/NONTRIPLET_DEG/NON_TRIPLET_DEG_down_tss.bed", header = FALSE) %>% 
  dplyr::select(V4)

pdf("../figures/MB03_6h_GEP_SIG_DEG_BINDINGPATTERNS_NON_TripletBound_overlapViolins.pdf", height = 5, width = 8)

nonTriplet_up %>% 
  merge(all_quant_pergene2, by = "V4") %>% 
  dplyr::select(-c(matches("RNA"))) %>%
  return_faltviolin_dmsoMs275()+
  ggtitle("MS275-UP, n=98/397; non-triplet")

nonTriplet_down %>% 
  merge(all_quant_pergene2, by = "V4") %>% 
  dplyr::select(-c(matches("RNA"))) %>%
  return_faltviolin_dmsoMs275()+
  ggtitle("MS275-DOWN, n=39/301; non-triplet")

dev.off()

## theme classic
pdf("../figures/MB03_6h_GEP_SIG_DEG_BINDINGPATTERNS_NON_TripletBound_overlapViolins_classic.pdf", height = 4.5, width = 6)

nonTriplet_up %>% 
  merge(all_quant_pergene2, by = "V4") %>% 
  dplyr::select(-c(matches("RNA"))) %>%
  return_faltviolin_dmsoMs275()+
  ggtitle("MS275-UP, n=98/397; non-triplet")+
  theme_classic(base_size = 20)+
  xlab("")+
  scale_x_discrete(expand = expand_scale(add = 0.35))+
  theme(axis.text = element_text(color = "black"),
        axis.ticks.length = unit(0.3, "cm"))

nonTriplet_down %>% 
  merge(all_quant_pergene2, by = "V4") %>% 
  dplyr::select(-c(matches("RNA"))) %>%
  return_faltviolin_dmsoMs275()+
  ggtitle("MS275-DOWN, n=39/301; non-triplet")+
  theme_classic(base_size = 20)+
  xlab("")+
  scale_x_discrete(expand = expand_scale(add = 0.35))+
  theme(axis.text = element_text(color = "black"),
        axis.ticks.length = unit(0.3, "cm"))

dev.off()


## overlapping violins of expression plots
return_gep_flatviolin <- function(gep_object, up_df, down_df){
  
  gep_mean <- gep_object %>% 
    dplyr::mutate(DMSO = c(DMSO_1+DMSO_2+DMSO_3)/3,
                  MS275 = c(MS275_1+MS275_2+MS275_3)/3) %>% 
    dplyr::select(GENE, DMSO, MS275) 
  
  up_exp <- up_df %>% 
    dplyr::select(V4) %>% 
    merge(gep_mean, by.x = "V4", by.y = "GENE") %>% 
    reshape2::melt() %>% 
    dplyr::mutate(group = "up")
  
  down_exp <- down_df %>% 
    dplyr::select(V4) %>% 
    merge(gep_mean, by.x = "V4", by.y = "GENE") %>% 
    reshape2::melt() %>% 
    dplyr::mutate(group = "down")
  
  # all in one
  all_exp <- rbind(up_exp, down_exp) %>% 
    dplyr::mutate(colr = paste0(variable, "_", group)) %>% 
    dplyr::mutate(colr = dplyr::case_when(colr == "DMSO_up" ~ "DMSO",
                                          colr == "DMSO_down" ~ "DMSO",
                                          colr == "MS275_down" ~ "MS275_down",
                                          colr == "MS275_up" ~ "MS275_up"))
  
  gep_cols <- c("blue", "red", "gray")
  names(gep_cols) <- c("MS275_down", "MS275_up", "DMSO")
  
  plt <- ggplot(data = all_exp, aes(y = value, x = group, fill = colr)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8)+
    geom_boxplot(width = .3, outlier.shape = NA, alpha = 0.8) +
    scale_color_manual(values = gep_cols)+
    scale_fill_manual(values = gep_cols)+
    theme_vt_ms()+
    ylab("Gene expression")+
    xlab("")+
    theme(legend.title = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank())+
    stat_compare_means(aes(group = variable), label = "p.format", size = 6, label.y = max(all_exp$value) + 0.5)
  
  return(plt)  
  
}

load("../data/HDMB03_GEP_6h_24h_48h_allReplicates_normalizedExpression.RData")

mb03_6h_up_all <- read.delim("~/Desktop/ServerView/MYC/results/affy-ge/DEG_FDR01_LOGFC1/ALL_UP_TSS.txt", header = FALSE)
mb03_6h_down_all <- read.delim("~/Desktop/ServerView/MYC/results/affy-ge/DEG_FDR01_LOGFC1/ALL_DOWN_TSS.txt", header = FALSE)

mb03_6h_up_triplet <- read.delim("~/Desktop/ServerView/MYC/results/affy-ge/DEG_FDR01_LOGFC1/TRIPLET_UP_TSS.txt", header = FALSE)
mb03_6h_down_triplet <- read.delim("~/Desktop/ServerView/MYC/results/affy-ge/DEG_FDR01_LOGFC1/TRIPLET_DOWN_TSS.txt", header = FALSE)


## print plots
pdf("../figures/HDMB03_6h_GEP_SIG_expressionPatterns_allTimePoints.pdf", height = 5, width = 8.5)

return_gep_flatviolin(gep_object = mb03_6h, up_df = mb03_6h_up_all, down_df = mb03_6h_down_all)+
  ggtitle("MB03-6h, all deg")

return_gep_flatviolin(gep_object = mb03_24h, up_df = mb03_6h_up_all, down_df = mb03_6h_down_all)+
  ggtitle("MB03-24h, all deg")

return_gep_flatviolin(gep_object = mb03_48h, up_df = mb03_6h_up_all, down_df = mb03_6h_down_all)+
  ggtitle("MB03-48h, all deg")

## tripelt only
return_gep_flatviolin(gep_object = mb03_6h, up_df = mb03_6h_up_triplet, down_df = mb03_6h_down_triplet)+
  ggtitle("MB03-6h, triplet deg")

return_gep_flatviolin(gep_object = mb03_24h, up_df = mb03_6h_up_triplet, down_df = mb03_6h_down_triplet)+
  ggtitle("MB03-24h, triplet deg")

return_gep_flatviolin(gep_object = mb03_48h, up_df = mb03_6h_up_triplet, down_df = mb03_6h_down_triplet)+
  ggtitle("MB03-48h, triplet deg")

dev.off()

## triplet bound DEG expression pattern before treatment (only DMSO)
mb03_6h_up_triplet_dmsoexp <- mb03_6h_up_triplet %>% 
  dplyr::select(V4) %>% 
  merge(mb03_6h, by.x = "V4", by.y = "GENE") %>% 
  dplyr::mutate(DMSO = (DMSO_1 + DMSO_2 + DMSO_3)/3,
                MS275 = (MS275_1 + MS275_2 + MS275_3)/3) %>% 
  dplyr::select(V4, DMSO) %>% 
  dplyr::mutate(group = "up") 

mb03_6h_updown_dmsoexp <- mb03_6h_down_triplet %>% 
  dplyr::select(V4) %>% 
  merge(mb03_6h, by.x = "V4", by.y = "GENE") %>% 
  dplyr::mutate(DMSO = (DMSO_1 + DMSO_2 + DMSO_3)/3,
                MS275 = (MS275_1 + MS275_2 + MS275_3)/3) %>% 
  dplyr::select(V4, DMSO) %>% 
  dplyr::mutate(group = "down") %>% 
  rbind(mb03_6h_up_triplet_dmsoexp) 
  

mb03_6h_updown_dmsoexp_violins <- ggplot(data = mb03_6h_updown_dmsoexp, aes(y = DMSO, x = group, fill = group)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8)+
  geom_boxplot(width = .2, outlier.shape = NA, alpha = 0.8) +
  scale_color_manual(values = c("down" = "blue", "up" = "red"))+
  scale_fill_manual(values = c("down" = "blue", "up" = "red"))+
  theme_vt_ms()+
  ylab("Gene expression")+
  xlab("")+
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  stat_compare_means(label = "p.format", size = 6, label.y = max(mb03_6h_updown_dmsoexp$DMSO) + 0.5)

mb03_6h_updown_dmsoexp_violins2 <- ggplot(data = mb03_6h_updown_dmsoexp, aes(y = DMSO, x = group, fill = group)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8)+
  geom_boxplot(width = .2, outlier.shape = NA, alpha = 0.8) +
  scale_color_manual(values = c("down" = "blue", "up" = "red"))+
  scale_fill_manual(values = c("down" = "blue", "up" = "red"))+
  theme_classic(base_size = 20)+
  ylab("Gene expression")+
  xlab("")+
  scale_x_discrete(expand = expand_scale(add = 0.2))+
  theme(axis.text = element_text(color = "black"),
        axis.ticks.length = unit(0.3, "cm"))+
  stat_compare_means(label = "p.format", size = 6, label.y = max(mb03_6h_updown_dmsoexp$DMSO) + 0.5)


save_plot("../figures/DEG_UP_DOWN_DMSO_ExpPatterns.pdf", mb03_6h_updown_dmsoexp_violins, base_aspect_ratio = 1.3, base_height = 5, base_width = 7)
save_plot("../figures/DEG_UP_DOWN_DMSO_ExpPatterns_classic.pdf", mb03_6h_updown_dmsoexp_violins2, base_aspect_ratio = 1.3, base_height = 4, base_width = 6)


### RNApolII metagene plot
plotProfile_plot_RNApolII_themeClassic <- function(melted_dat){
  
  # melted_dat <- up_rnapol_mlt
  
  ## custom labels
  lab_pos <- c(100, 400)
  lab_nm <- c("TSS", "TES")
  
  spline_int <- as.data.frame(spline(melted_dat$bin, melted_dat$value))
  
  my_plt <- ggplot(melted_dat, aes(bin, value, color = variable))+
    geom_line(size = 2)+
    scale_x_continuous(breaks = lab_pos, labels = lab_nm, expand = expand_scale(mult = c(0, 0)))+
    theme_vt_ms()+
    xlab("")+
    ylab("input normalized RNApolII signal (Log2)")+
    #scale_color_manual(values = c("RNApolII_untreated" = "darkgray", "RNApolII_treated" = "black"))+
    theme(legend.title = element_blank(),
          legend.position = c(0.85, 0.9),
          panel.grid = element_blank())
  return(my_plt)
  
}

TRIPLET_RNAPOL <- "~/Desktop/ServerView/MYC/results/affy-ge/DEG_FDR01_LOGFC1"

### only before treatment (red,blue)
dmso_up_down_RNApolII = paste0(TRIPLET_RNAPOL, "/HDMB03_TRIPLET_DEG_RNApolII.data")

return_plotProfile_formated_RNApolII_dmso <- function(file_path){
  
  # file_path = dmso_up_down_RNApolII
  dat <- read.delim(file_path, header = TRUE, skip = 1)
  
  dat_mlt <- dat %>%  
    dplyr::select(-c(bins)) %>% 
    tibble::column_to_rownames(var = "X") %>% t() %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "bin") %>% 
    dplyr::filter(bin != "X") %>% 
    dplyr::mutate(bin = gsub("X", "", bin)) %>% 
    reshape2::melt() %>% 
    dplyr::mutate(bin = as.numeric(bin)) 
  
  return(dat_mlt)
  
}

dmso_up_down_RNApolII_mlt <- return_plotProfile_formated_RNApolII_dmso(file_path = dmso_up_down_RNApolII)

pdf("../figures/HDMB03_6h_GEP_SIG_DEG_Triplet_RNApolII_metagene_beforeTreatment.pdf", height = 5.5, width = 8)

plotProfile_plot_RNApolII_themeClassic(melted_dat = dmso_up_down_RNApolII_mlt)+
  scale_color_manual(values = c("up" = "red", "down" = "blue"))

dev.off()

### triplet RNAPolII change before and after treatment
triplet_down_rnapol <- paste0(TRIPLET_RNAPOL, "/TRIPLET_DOWN_wholegene_RNApolII.data")
triplet_up_rnapol <- paste0(TRIPLET_RNAPOL, "/TRIPLET_UP_wholegene_RNApolII.data")

return_plotProfile_formated_RNApolII_Un_Tr <- function(file_path){
  
  # file_path = dmso_up_down_RNApolII
  dat <- read.delim(file_path, header = TRUE, skip = 1)
  
  dat_mlt <- dat %>%  
    dplyr::select(-c(X)) %>% 
    tibble::column_to_rownames(var = "bins") %>% t() %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "bin") %>% 
    dplyr::filter(bin != "X") %>% 
    dplyr::mutate(bin = gsub("X", "", bin)) %>% 
    reshape2::melt() %>% 
    dplyr::mutate(bin = as.numeric(bin)) 
  
  return(dat_mlt)
  
}

triplet_down_rnapol_mlt <- return_plotProfile_formated_RNApolII_Un_Tr(file_path = triplet_down_rnapol)
triplet_up_rnapol_mlt <- return_plotProfile_formated_RNApolII_Un_Tr(file_path = triplet_up_rnapol)

pdf("../figures/HDMB03_6h_GEP_SIG_DEG_Triplet_RNApolII_metagene_DMSO_MS275_classic.pdf", height = 4.5, width = 6)

plotProfile_plot_RNApolII_themeClassic(melted_dat = triplet_down_rnapol_mlt)+
  scale_color_manual(values = c("RNApolII_untreated" = "darkgray", "RNApolII_treated" = "black"))+
  theme_classic(base_size = 18)+
  theme(axis.text = element_text(color = "black"),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_blank(),
        legend.position = c(0.7, 0.9))+
  ggtitle("Triplet-down")

plotProfile_plot_RNApolII_themeClassic(melted_dat = triplet_up_rnapol_mlt)+
  scale_color_manual(values = c("RNApolII_untreated" = "darkgray", "RNApolII_treated" = "black"))+
  theme_classic(base_size = 18)+
  theme(axis.text = element_text(color = "black"),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_blank(),
        legend.position = c(0.7, 0.9))+
  ggtitle("Triplet-up")

dev.off()

##
pdf("../figures/HDMB03_6h_GEP_SIG_DEG_Triplet_RNApolII_metagene_beforeTreatment_classic.pdf", height = 4.5, width = 6)

plotProfile_plot_RNApolII_themeClassic(melted_dat = dmso_up_down_RNApolII_mlt)+
  scale_color_manual(values = c("up" = "red", "down" = "blue"))+
  theme_classic(base_size = 18)+
  scale_y_continuous(breaks = c(0, 1, 2, 3), labels = c("0", "1", "2", "3"))+
  theme(axis.text = element_text(color = "black"),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_blank(),
        legend.position = c(0.7, 0.9))

dev.off()

### before treatment, protein binding patterns on DEG

all_quant2_dmso = all_quant_pergene2 %>% 
  dplyr::select(matches("V4|_untreated")) %>% 
  dplyr::select(-c(RNApolII_untreated)) 

colnames(all_quant2_dmso) = colnames(all_quant2_dmso) %>% {gsub("_untreated", "", .)}
all_quant2_dmso = all_quant2_dmso %>% dplyr::select(V4, MYC, HDAC2, H3K27ac)

cl_triplet_down_clQuant = cl_triplet_down %>% 
  dplyr::select(GENE) %>% 
  merge(all_quant2_dmso, by.x = "GENE", by.y = "V4") %>% 
  #dplyr::select(-c(RNApolII)) %>% 
  reshape2::melt() %>% 
  dplyr::mutate(DEG = "DOWN") 

cl_triplet_up_clQuant = cl_triplet_up %>% 
  dplyr::select(GENE) %>% 
  merge(all_quant2_dmso, by.x = "GENE", by.y = "V4") %>% 
  #dplyr::select(-c(RNApolII)) %>% 
  reshape2::melt() %>% 
  dplyr::mutate(DEG = "UP")

cl_upDown_clQuant = rbind(cl_triplet_down_clQuant, cl_triplet_up_clQuant) %>% 
  dplyr::mutate(DEG = dplyr::case_when(DEG == "DOWN" ~ "down",
                                       DEG == "UP" ~ "up"))

pdf("../figures/HDMB03_6h_GEP_SIG_DEG_BINDINGPATTERNS_TripletBound_BeforeTreatment_overlapViolins.pdf", height = 5, width = 8)

ggplot(data = cl_upDown_clQuant, aes(y = value, x = variable, fill = DEG)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .7)+
  #geom_point(aes(y = value, color = DEG), position = position_jitter(width = .05), size = 1, alpha = 0.8) +
  geom_boxplot(width = .3, outlier.shape = NA) +
  #guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_fill_manual(values = c("blue", "red"))+
  scale_color_manual(values = c("blue", "red"))+
  theme_vt_ms()+
  ylab("input normalized IP signal (log2)")+
  stat_compare_means(aes(group = DEG), label = "p.format", size = 6, label.y = max(cl_upDown_clQuant$value) + 0.5)+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())+
  # legend.position = c(0.1, 0.65),
  # legend.background = element_rect(fill=NA,
  #                                  size=0.5, linetype="solid", 
  #                                  colour ="black"))+
  ggtitle("HD-MB03 Triplet DEG - DMSO\n(Down-261; Up-297)")

ggplot(data = cl_upDown_clQuant, aes(y = value, x = variable, fill = DEG)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .7)+
  #geom_point(aes(y = value, color = DEG), position = position_jitter(width = .05), size = 1, alpha = 0.8) +
  geom_boxplot(width = .3, outlier.shape = NA) +
  #guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_fill_manual(values = c("blue", "red"))+
  scale_color_manual(values = c("blue", "red"))+
  theme_vt_ms()+
  ylab("input normalized IP signal (log2)")+
  stat_compare_means(aes(group = DEG), label = "p.format", size = 6, label.y = max(cl_upDown_clQuant$value) + 0.5)+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())+
        # legend.position = c(0.1, 0.65),
        # legend.background = element_rect(fill=NA,
        #                                  size=0.5, linetype="solid", 
        #                                  colour ="black"))+
  ggtitle("HD-MB03 Triplet DEG - DMSO (Down-261; Up-297)")

dev.off()

###

pdf("../figures/HDMB03_6h_GEP_SIG_DEG_BINDINGPATTERNS_TripletBound_BeforeTreatment_overlapViolins_classic.pdf", height = 4.5, width = 7)

ggplot(data = cl_upDown_clQuant, aes(y = value, x = variable, fill = DEG)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .7)+
  #geom_point(aes(y = value, color = DEG), position = position_jitter(width = .05), size = 1, alpha = 0.8) +
  geom_boxplot(width = .3, outlier.shape = NA) +
  #guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_fill_manual(values = c("blue", "red"))+
  scale_color_manual(values = c("blue", "red"))+
  theme_classic(base_size = 20)+
  ylab("input normalized IP signal (log2)")+
  stat_compare_means(aes(group = DEG), label = "p.format", size = 6, label.y = max(cl_upDown_clQuant$value) + 0.5)+
  theme_classic(base_size = 20)+
  xlab("")+
  scale_x_discrete(expand = expand_scale(add = 0.35))+
  theme(axis.text = element_text(color = "black"),
        axis.ticks.length = unit(0.3, "cm"))
  #ggtitle("HD-MB03 Triplet DEG - DMSO\n(Down-261; Up-297)")

dev.off()

### global changes in protein binding after treatment, ecdf
## all proteins are quantified at TSS/TES with `deeptools multiBigwigSummary` function

d_RNAPolII <- readRDS("../data/HDMB03_DMSO_MS275_RNApolII_TES.rds")
d_dmso_gt0 = d_RNAPolII %>% 
  dplyr::filter(RNApolII_untreated > 0 & RNApolII_treated > 0)

d_mlt_all <- d_RNAPolII %>% dplyr::mutate(REGION = paste(chr, start, end, sep = "_")) %>% 
  dplyr::select(-c(chr:end)) %>% 
  reshape2::melt()

d_mlt = d_dmso_gt0 %>% dplyr::mutate(REGION = paste(chr, start, end, sep = "_")) %>% 
  dplyr::select(-c(chr:end)) %>% 
  reshape2::melt() 

## all proteins, melted before saving from original script
all_proteins <- readRDS("../data/HDMB03_DMSO_MS275_allProteins_TSS.rds")

myc = all_proteins %>% dplyr::filter(grepl("MYC", variable))
hdac2 = all_proteins %>% dplyr::filter(grepl("HDAC2", variable))
h3k27ac = all_proteins %>% dplyr::filter(grepl("H3K27ac", variable))

pdf("../figures/AllProteins_AllProteinCodingGenes_DMSO_MS275_ecdf.pdf", height = 7, width = 6.5)

## TES RNApolII
ggplot(d_mlt_all, aes(x = value, color = variable))+
  stat_ecdf(size = 1.5)+
  scale_color_manual(values = c("black","darkgray"))+
  scale_y_continuous(labels = scales::percent)+
  theme_vt_ms()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.7, 0.15),
        legend.background = element_rect(fill= "white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"))+
  ylab("% of protein coding genes")+
  xlab("input normalized TES RNApolII signal (Log2)")

## TSS MYC
ggplot(myc, aes(x = value, color = variable))+
  stat_ecdf(size = 1.5)+
  scale_color_manual(values = c("#78CA20","darkgray"))+
  scale_y_continuous(labels = scales::percent)+
  theme_vt_ms()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.7, 0.15),
        legend.background = element_rect(fill= "white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"))+
  ylab("% of protein coding genes")+
  xlab("input normalized MYC signal (Log2)")

## TSS HDAC2
ggplot(hdac2, aes(x = value, color = variable))+
  stat_ecdf(size = 1.5)+
  scale_color_manual(values = c("#AD0AFD","darkgray"))+
  scale_y_continuous(labels = scales::percent)+
  theme_vt_ms()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.7, 0.15),
        legend.background = element_rect(fill= "white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"))+
  ylab("% of protein coding genes")+
  xlab("input normalized HDAC2 signal (Log2)")

## TSS H3K27ac
ggplot(h3k27ac, aes(x = value, color = variable))+
  stat_ecdf(size = 1.5)+
  scale_color_manual(values = c("#F22C1E","darkgray"))+
  scale_y_continuous(labels = scales::percent)+
  theme_vt_ms()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.7, 0.15),
        legend.background = element_rect(fill= "white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"))+
  ylab("% of protein coding genes")+
  xlab("input normalized H3K27ac signal (Log2)")

dev.off()

##
pdf("../figures/AllProteins_AllProteinCodingGenes_DMSO_MS275_ecdf_classic.pdf", height = 6, width = 5.6)

## TES RNApolII
ggplot(d_mlt_all, aes(x = value, color = variable))+
  stat_ecdf(size = 1.5)+
  scale_color_manual(values = c("black","darkgray"))+
  scale_y_continuous(labels = scales::percent)+
  theme_vt_ms()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.7, 0.15))+
  ylab("% of protein coding genes")+
  xlab("input normalized TES RNApolII signal (Log2)")

## TSS MYC
ggplot(myc, aes(x = value, color = variable))+
  stat_ecdf(size = 1.5)+
  scale_color_manual(values = c("#78CA20","darkgray"))+
  scale_y_continuous(labels = scales::percent)+
  theme_vt_ms()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.7, 0.15))+
  ylab("% of protein coding genes")+
  xlab("input normalized MYC signal (Log2)")

## TSS HDAC2
ggplot(hdac2, aes(x = value, color = variable))+
  stat_ecdf(size = 1.5)+
  scale_color_manual(values = c("#AD0AFD","darkgray"))+
  scale_y_continuous(labels = scales::percent)+
  theme_vt_ms()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.7, 0.15))+
  ylab("% of protein coding genes")+
  xlab("input normalized HDAC2 signal (Log2)")

## TSS H3K27ac
ggplot(h3k27ac, aes(x = value, color = variable))+
  stat_ecdf(size = 1.5)+
  scale_color_manual(values = c("#F22C1E","darkgray"))+
  scale_y_continuous(labels = scales::percent)+
  theme_vt_ms()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.7, 0.15))+
  ylab("% of protein coding genes")+
  xlab("input normalized H3K27ac signal (Log2)")

dev.off()

## determine p.values

pdf("../figures/MYC_HDAC2_H3K27ac_tss_change_global_boxpvals.pdf", height = 6, width = 4.5)

myc %>% 
  ggpubr::ggboxplot(x = "variable", y = "value")+
  stat_compare_means(size = 6)

hdac2 %>% 
  ggpubr::ggboxplot(x = "variable", y = "value")+
  stat_compare_means(size = 6)

h3k27ac %>% 
  ggpubr::ggboxplot(x = "variable", y = "value")+
  stat_compare_means(size = 6)

dev.off()



################################################
################################################
## GSEA plots for 3 different cell-lines
library(fgsea)

EXP.INPUT <- "~/Desktop/ServerView/MYC/results/affy-ge/all_timepoints_deg"

mb03_exp <- read.delim(paste0(EXP.INPUT, "/ALL_DEG_6h.txt"), header = TRUE)
d458_exp <- read.delim(paste0(EXP.INPUT, "/ALL_DEG_D458_6h.txt"), header = TRUE)
med8a_exp <- read.delim(paste0(EXP.INPUT, "/ALL_DEG_MED8A_6h.txt"), header = TRUE)

hallmark_v2 <- read.delim("../data/HALLMARK_MYC_TARGETS_V2.txt", header = FALSE)
hallmark_v2_set <- list()
hallmark_v2_set[['HALLMARK_MYC_TARGETS_V2']] <- hallmark_v2$V1 %>% as.character() %>% unique()

plotEnrichment_custom <- function(pathway, stats,
                           gseaParam=1,
                           ticksSize=0.2) {
  
  rnk <- rank(-stats)
  ord <- order(rnk)
  
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))
  
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                          returnAllExtremes = TRUE)
  
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))
  
  diff <- (max(tops) - min(bottoms)) / 8
  
  # Getting rid of NOTEs
  x=y=NULL
  g <- ggplot(toPlot, aes(x=x, y=y)) +
    geom_point(color="green", size=0.3) +
    geom_hline(yintercept=max(tops), colour="red", linetype="dashed", size = 0.9) +
    geom_hline(yintercept=min(bottoms), colour="red", linetype="dashed", size = 0.9) +
    geom_hline(yintercept=0, colour="black", size = 0.9) +
    geom_line(color="green", size = 0.9) + theme_bw() +
    geom_segment(data=data.frame(x=pathway),
                 mapping=aes(x=x, y=-diff/2,
                             xend=x, yend=diff/2),
                 size=0.6) +
    
    theme(panel.border=element_blank(),
          panel.grid.minor=element_blank()) +
    
    labs(x="rank", y="enrichment score")
  g
}

return_gsea_plot <- function(deg_df, gene_set = hallmark_v2_set){
  
  my_ranks <- deg_df %>% 
    dplyr::select(GENE, logFC) %>% 
    dplyr::arrange(desc(logFC)) %>% 
    tibble::deframe()  
  
  fgseaRes <- fgsea(pathways=gene_set, stats=my_ranks, nperm=1000)
  
  plotEnrichment_custom(gene_set[["HALLMARK_MYC_TARGETS_V2"]],
                 my_ranks) + labs(title="HALLMARK_MYC_TARGETS_V2")+
    theme_classic(base_size = 18)+
    theme(axis.ticks.length = unit(0.3, "cm"),
          axis.text = element_text(color = "black"))
}

## MB03 - NES: -1.18
## D458 - NES: -1.48
## MED8A - NES: -1.83

pdf("../figures/HALLMARK_MYC_TARGETS_V2_gseaplot.pdf", height = 4, width = 5.5)

return_gsea_plot(deg_df = mb03_exp)+ggtitle("HDMB03")
return_gsea_plot(deg_df = d458_exp)+ggtitle("D458")
return_gsea_plot(deg_df = med8a_exp)+ggtitle("MED8A")

dev.off()





