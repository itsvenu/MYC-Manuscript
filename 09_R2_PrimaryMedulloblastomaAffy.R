## source: /b06x-lsdf/infra5-lsdf/mbHDAC2/scripts/data_290917/Figures/R2_MYC/R2_MYC_vs_TARGETS_MB.tmp.R

## R2 genomics MB Affy cohort
## With increased MYC expression, do we have overall increased 
## `HALLMARK_MYC_TARGET_V2` gene expression in any of the subgroups?

library(tidyverse)

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
          panel.border = element_rect(size=1.6, color="black", fill=NA),
          axis.text = element_text(color="black"),
          axis.title = element_text(color="black"),
          legend.text = element_text(color="black"),
          legend.title = element_text(color="black")
    )
}

v2_dat = readRDS("../data/R2_Medulloblastoma_Hallmark_MYC_targetGenes_AffyExp.rds")
all_in = readRDS("../data/R2_Medulloblastoma_Hallmark_MYC_targetGenes_AffyExp_V1V2.rds")

return_loess_plot_v02 <- function(subgroup_name, hex_color){
  
  my_pal = hex_color
  names(my_pal) = subgroup_name
  
  my_myc_ordered = v2_dat %>% 
    dplyr::filter(subgroup == get("subgroup_name")) %>% 
    dplyr::arrange(MYC) %>% 
    dplyr::select(-c(MYC, subgroup)) %>% 
    tibble::column_to_rownames(var = "PID")
  
  ## total number of samples
  my_total_pid = dim(my_myc_ordered)[1]
  
  my_myc_ordered[] <- lapply(my_myc_ordered, function(x) {
    if(is.factor(x)) as.numeric(as.character(x)) else x
  })
  
  my_targetMean_MYC = rowMeans(my_myc_ordered) %>% 
    as.data.frame() %>% 
    dplyr::rename(TARGETS_MEAN = ".") %>% 
    tibble::rownames_to_column(var = "PID") %>% 
    merge(all_in %>% 
            dplyr::select(PID, MYC))
  
  my_targetMean_MYC$MYC = as.numeric(as.character(my_targetMean_MYC$MYC))
  
  my_targetMean_MYC = my_targetMean_MYC %>% 
    dplyr::arrange(MYC) %>% 
    tibble::rowid_to_column(var = "ID") %>% 
    dplyr::mutate(SUBGROUP = subgroup_name)
  
  my_loess_plot = ggplot(my_targetMean_MYC, aes(MYC, TARGETS_MEAN, color = SUBGROUP))+
    geom_smooth(method = "loess", se = FALSE, size = 1.7)+
    theme_vt()+
    scale_color_manual(values = my_pal)+
    ylab("HALLMARK_MYC_TARGETS_V2 (mean)")+
    xlab("Ranked MYC expression")+
    theme(legend.position = "none",
          panel.grid = element_blank())+
    ggtitle(paste0(subgroup_name, ", n=", my_total_pid))
  
  ## corr plot
  # tt = my_targetMean_MYC %>% 
  #   dplyr::filter(PID != "sjmb001")
  
  my_corr_plot = ggscatter(my_targetMean_MYC, x = "MYC", y = "TARGETS_MEAN", color = "SUBGROUP", palette = my_pal, size = 3)+
    stat_cor(method = "pearson", size = 6)+
    geom_smooth(method = "loess", se = FALSE, size = 1.7, color = "black")+
    theme_vt()+
    ylab("HALLMARK_MYC_TARGETS_V2 (mean)")+
    xlab("Ranked MYC expression")+
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank())+
    ggtitle(paste0(subgroup_name, ", n=", my_total_pid))
  
  pdf(paste0("../figures/HALLMARK_MYC_TARGETS_V2_vs_MYC_", subgroup_name, "_ExpPlot.pdf"), height = 5, width = 6.5)
  print(my_loess_plot)
  print(my_corr_plot)
  dev.off()
  
  return(my_targetMean_MYC)
  
  message(" ## DONE: ", subgroup_name)
}

## plot
group3_exp <- return_loess_plot_v02(subgroup_name = "group_3", hex_color = "#FFC000FF")
group4_exp <- return_loess_plot_v02(subgroup_name = "group_4", hex_color = "#4F6228FF")
ssh_exp <- return_loess_plot_v02(subgroup_name = "shh", hex_color = "#C00000FF")
wnt_exp <- return_loess_plot_v02(subgroup_name = "wnt", hex_color = "#004586FF")

all_subgroups <- rbind(group3_exp, group4_exp, ssh_exp, wnt_exp) %>% 
  dplyr::mutate(sub_group = case_when(SUBGROUP == "group_3" ~ "Group-3 (n=56)",
                                      SUBGROUP == "group_4" ~ "Group-4 (n=91)",
                                      SUBGROUP == "shh" ~ "SHH (n=59)",
                                      SUBGROUP == "wnt" ~ "WNT (n=17)")) 

subgroupc_colors <- c("#FFC000FF", "#4F6228FF", "#C00000FF", "#004586FF")
names(subgroupc_colors) <- c("group_3", "group_4", "shh", "wnt")

my_plt <- ggplot(all_subgroups, aes(MYC, TARGETS_MEAN, color = SUBGROUP))+
  geom_point()+
  theme_vt()+
  scale_color_manual(values = subgroupc_colors)+
  stat_cor(size = 6, label.y = 8.5, color = "black")+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour="black", fill=NA, size = 1.3),
        legend.position = "none")+
  facet_wrap(.~sub_group, scales = "free_x", ncol = 2)+
  xlab("MYC expression")+
  ylab("Hallmark MYC target gene expression (mean)")

all_subgroups_classic_plt <- ggplot(all_subgroups, aes(MYC, TARGETS_MEAN, color = SUBGROUP))+
  geom_point()+
  geom_smooth(method = "loess", se = FALSE, color = "black", size = 0.9)+
  scale_color_manual(values = subgroupc_colors)+
  stat_cor(size = 6, label.y = 8.5, color = "black")+
  facet_wrap(.~sub_group, scales = "free_x", ncol = 2)+
  theme_classic(base_size = 18)+
  theme(strip.background = element_rect(colour=NA, fill=NA),
        legend.position = "none",
        axis.ticks.length = unit(0.3, "cm"),
        axis.text = element_text(color = "black"),
        strip.text.x = element_text(size = 16))+
  xlab("MYC expression")+
  ylab("Hallmark MYC target genes expression (mean)")


all_subgroups_classic_plt2 <- ggplot(all_subgroups, aes(MYC, TARGETS_MEAN, color = SUBGROUP))+
  geom_point()+
  geom_smooth(method = "loess", se = FALSE, color = "black", size = 0.9)+
  scale_color_manual(values = subgroupc_colors)+
  ggpubr::stat_cor(size = 6, label.y = 8.5, color = "black")+
  facet_wrap(.~sub_group, scales = "free_x", ncol = 4)+
  theme_classic(base_size = 18)+
  theme(strip.background = element_rect(colour=NA, fill=NA),
        legend.position = "none",
        axis.ticks.length = unit(0.3, "cm"),
        axis.text = element_text(color = "black"),
        strip.text.x = element_text(size = 16),
        panel.spacing = unit(3, "lines"))+
  xlab("MYC expression")+
  ylab("Hallmark MYC target genes\nexpression (mean)")


save_plot("../figures/HALLMARK_MYC_TARGETS_V2_vs_MYC_allSubgroups_classic.pdf", all_subgroups_classic_plt, base_aspect_ratio = 1.3, base_height = 7, base_width = 7.5)

save_plot("../figures/HALLMARK_MYC_TARGETS_V2_vs_MYC_allSubgroups_classic2.pdf", all_subgroups_classic_plt2, base_aspect_ratio = 1.3, base_height = 3.5, base_width = 14)


g3_plot <- ggplot(all_subgroups %>% 
         dplyr::filter(SUBGROUP == "group_3"), aes(MYC, TARGETS_MEAN, color = SUBGROUP))+
  geom_point()+
  geom_smooth(method = "loess", se = FALSE, color = "black", size = 0.9)+
  scale_color_manual(values = subgroupc_colors)+
  stat_cor(size = 6, label.y = 8.5, color = "black")+
  facet_wrap(.~sub_group, scales = "free_x", ncol = 2)+
  theme_classic(base_size = 18)+
  theme(strip.background = element_rect(colour=NA, fill=NA),
        legend.position = "none",
        axis.ticks.length = unit(0.3, "cm"),
        axis.text = element_text(color = "black"),
        strip.text.x = element_text(size = 16))+
  xlab("")+
  ylab("")



g4_plot <- ggplot(all_subgroups %>% 
                    dplyr::filter(SUBGROUP == "group_4"), aes(MYC, TARGETS_MEAN, color = SUBGROUP))+
  geom_point()+
  geom_smooth(method = "loess", se = FALSE, color = "black", size = 0.9)+
  scale_color_manual(values = subgroupc_colors)+
  stat_cor(size = 6, label.y = 8.5, color = "black")+
  facet_wrap(.~sub_group, scales = "free_x", ncol = 2)+
  theme_classic(base_size = 18)+
  theme(strip.background = element_rect(colour=NA, fill=NA),
        legend.position = "none",
        axis.ticks.length = unit(0.3, "cm"),
        axis.text = element_text(color = "black"),
        strip.text.x = element_text(size = 16))+
  xlab("")+
  ylab("")

shh_plot <- ggplot(all_subgroups %>% 
                    dplyr::filter(SUBGROUP == "shh"), aes(MYC, TARGETS_MEAN, color = SUBGROUP))+
  geom_point()+
  geom_smooth(method = "loess", se = FALSE, color = "black", size = 0.9)+
  scale_color_manual(values = subgroupc_colors)+
  stat_cor(size = 6, label.y = 8.5, color = "black")+
  facet_wrap(.~sub_group, scales = "free_x", ncol = 2)+
  theme_classic(base_size = 18)+
  theme(strip.background = element_rect(colour=NA, fill=NA),
        legend.position = "none",
        axis.ticks.length = unit(0.3, "cm"),
        axis.text = element_text(color = "black"),
        strip.text.x = element_text(size = 16))+
  xlab("")+
  ylab("")

wnt_plot <- ggplot(all_subgroups %>% 
                    dplyr::filter(SUBGROUP == "wnt"), aes(MYC, TARGETS_MEAN, color = SUBGROUP))+
  geom_point()+
  geom_smooth(method = "loess", se = FALSE, color = "black", size = 0.9)+
  scale_color_manual(values = subgroupc_colors)+
  stat_cor(size = 6, label.y = 8.5, color = "black")+
  facet_wrap(.~sub_group, scales = "free_x", ncol = 2)+
  theme_classic(base_size = 18)+
  theme(strip.background = element_rect(colour=NA, fill=NA),
        legend.position = "none",
        axis.ticks.length = unit(0.3, "cm"),
        axis.text = element_text(color = "black"),
        strip.text.x = element_text(size = 16))+
  xlab("")+
  ylab("")

gridExtra::grid.arrange(g3_plot, g4_plot, shh_plot, wnt_plot, ncol = 4)


# save_plot("../figures/HALLMARK_MYC_TARGETS_V2_vs_MYC_allSubgroups.pdf", my_plt, base_aspect_ratio = 1.3, base_height = 8, base_width = 8.5)

all_subgroups2 <- all_subgroups
all_subgroups2$sub_group <- factor(all_subgroups2$sub_group, levels = c("WNT (n=17)", "SHH (n=59)",
                                                                        "Group-3 (n=56)", "Group-4 (n=91)"))

all_subgroups_classic_corline_plt <- ggscatter(all_subgroups2, x = "MYC", y = "TARGETS_MEAN", color = "SUBGROUP",
          add = "reg.line", conf.int = FALSE,
          cor.coef = TRUE,
          cor.coeff.args = list(method = "pearson", label.sep = "\n"))+
  scale_color_manual(values = subgroupc_colors)+
  facet_wrap(.~sub_group, scales = "free_x", ncol = 4)+
  theme_classic(base_size = 18)+
  theme(strip.background = element_rect(colour=NA, fill=NA),
        legend.position = "none",
        axis.ticks.length = unit(0.3, "cm"),
        axis.text = element_text(color = "black"),
        strip.text.x = element_text(size = 16),
        panel.spacing = unit(3, "lines"))+
  xlab("MYC expression")+
  ylab("Hallmark MYC target genes\nexpression (mean)")

cowplot::save_plot("../figures/HALLMARK_MYC_TARGETS_V2_vs_MYC_allSubgroups_classic3_corline.pdf", all_subgroups_classic_corline_plt, base_aspect_ratio = 1.3, base_height = 3.5, base_width = 14)

  

