## source: /b06x-lsdf/infra5-lsdf/mbHDAC2/scripts/affy-ge/DEG_DMEA/E_BOX/ebox_metageneplot.R

## ranking eboxes
## source: /b06x-lsdf/infra5-lsdf/mbHDAC2/scripts/data_30052018/downstreamAnalysis/Untreated_treated_analysis/analysis_12102018/01_peak_filtering/E_box_ranking/rank_eboxes.py
## source: /b06x-lsdf/infra5-lsdf/mbHDAC2/scripts/data_30052018/downstreamAnalysis/Untreated_treated_analysis/analysis_12102018/01_peak_filtering/E_box_ranking/rankEboxes.R

## Ebox analysis on public Myc target genes from MsigDB
## source: /b06x-lsdf/infra5-lsdf/mbHDAC2/scripts/affy-ge/all_timepoints_deg/PUBLIC_MYC_TARGETS/msigDB_MYC_targets.R

## homer annotatePeaks function ran in `hist` mode
## using custom E-box motif files (prepared with `seq2profile.pl` function)

setwd("~/Desktop/ServerView/MYC/scripts/MycHdaci_code/scripts")

library(tidyverse)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

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

## read and format `homer annotatePeaks -hist` results
readHomer_metagene <- function(file_path){
  
  d = read.delim(file_path, header = TRUE)
  
  cl_names = colnames(d)
  cl_names[1] = "bin"
  
  colnames(d) = cl_names
  
  ## extract bin, *.total.sites* columns
  
  d = d %>% dplyr::select(matches("bin|.total.sites"))
  
  colnames(d) = d %>% colnames() %>% {gsub("\\.total.sites", "", .)}
  
  return(d)
  
}

## plot motif profile plots around TSS
ebox_metageneplot2 <- function(homer_df){
  
  ebox_meta = homer_df
  
  ebox_meta$bin = as.factor(ebox_meta$bin)
  
  ebox_meta_mlt = ebox_meta %>% reshape2::melt()
  
  #as.numeric(levels(f))[f]
  
  ebox_meta_mlt$bin = as.numeric(levels(ebox_meta_mlt$bin))[ebox_meta_mlt$bin]
  
  my_p = ggplot(ebox_meta_mlt, aes(bin, value, color = variable))+
    #geom_line(size = 1.5, alpha = 0.4)+
    geom_smooth(method = "loess", se = FALSE, size = 1.5)+
    facet_wrap(. ~ variable, ncol = 4, scales = "free_y")+
    theme_bw(base_size = 14)+
    theme(legend.position = "none")+
    ylab("motifs per bp per promoter")+
    xlab("")
  
  return(my_p)
}

## plot up/down ebox dist in one plot
return_ebox_overlapplot3 <- function(up_df, down_df, ebox_id){
  
  down_up_df = down_df %>% dplyr::select(bin, get("ebox_id")) %>% 
    dplyr::rename(down = get("ebox_id")) %>% 
    merge(up_df %>% dplyr::select(bin, get("ebox_id")) %>% 
            dplyr::rename(up = get("ebox_id")), by = "bin")
  
  down_up_df$bin = as.factor(down_up_df$bin)
  down_up_df_mlt = down_up_df %>% reshape2::melt()
  
  down_up_df_mlt$bin = as.numeric(levels(down_up_df_mlt$bin))[down_up_df_mlt$bin]
  
  my_loess_spd = ggplot(down_up_df_mlt, aes(x = bin, y = value, color = variable))+
    geom_smooth(method = "loess", size = 1.4, se  =FALSE)+
    scale_color_manual(values = c("down" = "blue", "up" = "red"))+
    theme_vt()+
    ggtitle(paste0(ebox_id, ", loess fit, spd"))+
    scale_x_continuous(breaks = c(-2000, -1000, 0, 1000, 2000), labels = c("-2000", "-1000", "TSS", "1000", "2000"))+
    theme(legend.title = element_blank(),
          panel.grid = element_blank())+
    ylab("motifs per bp per promoter")+
    xlab("distance from TSS")
  
  all_plt <- my_loess_spd
  
  return(all_plt)
  
}

return_ebox_overlapplot4 <- function(up_df, down_df, ebox_id){
  
  down_up_df = down_df %>% dplyr::select(bin, get("ebox_id")) %>% 
    dplyr::rename(down = get("ebox_id")) %>% 
    merge(up_df %>% dplyr::select(bin, get("ebox_id")) %>% 
            dplyr::rename(up = get("ebox_id")), by = "bin")
  
  down_up_df$bin = as.factor(down_up_df$bin)
  down_up_df_mlt = down_up_df %>% reshape2::melt()
  
  down_up_df_mlt$bin = as.numeric(levels(down_up_df_mlt$bin))[down_up_df_mlt$bin]
  
  my_loess_spd = ggplot(down_up_df_mlt, aes(x = bin, y = value, color = variable))+
    geom_smooth(method = "loess", size = 1.4, se  =FALSE)+
    scale_color_manual(values = c("down" = "blue", "up" = "red"))+
    theme_vt_ms()+
    ggtitle(paste0(ebox_id, ", loess fit, spd"))+
    scale_x_continuous(breaks = c(-2000, -1000, 0, 1000, 2000), labels = c("-2000", "-1000", "TSS", "1000", "2000"))+
    theme(legend.title = element_blank(),
          panel.grid = element_blank())+
    ylab("motifs per bp per promoter")+
    xlab("distance from TSS")
  
  all_plt <- my_loess_spd
  
  return(all_plt)
  
}

TRIPLET_UP <- "../data/TRIPLET_UP_TSS_10bp_eBoxDesnities.txt"
TRIPLET_DOWN <- "../data/TRIPLET_DOWN_TSS_10bp_eBoxDesnities.txt"

genomewide_ebox_nums <- readRDS("../data/Genomewide_eboxHits_perMillions.rds")

TRIPLET_UP_dist <- readHomer_metagene(file_path = TRIPLET_UP)
colnames(TRIPLET_UP_dist) <- gsub("EBOX_", "", colnames(TRIPLET_UP_dist))

TRIPLET_UP_dist2 <- TRIPLET_UP_dist %>% tibble::column_to_rownames(var = "bin")

TRIPLET_UP_dist2_norm <- TRIPLET_UP_dist2/genomewide_ebox_nums$pm_hits[match(names(TRIPLET_UP_dist2), genomewide_ebox_nums$ebox)][col(TRIPLET_UP_dist2)]
TRIPLET_UP_dist2_norm <- TRIPLET_UP_dist2_norm %>% tibble::rownames_to_column(var = "bin")
TRIPLET_UP_dist2_norm$bin <- as.numeric(TRIPLET_UP_dist2_norm$bin)

## down
TRIPLET_DOWN_dist <- readHomer_metagene(file_path = TRIPLET_DOWN)
colnames(TRIPLET_DOWN_dist) <- gsub("EBOX_", "", colnames(TRIPLET_DOWN_dist))

TRIPLET_DOWN_dist2 <- TRIPLET_DOWN_dist %>% tibble::column_to_rownames(var = "bin")

TRIPLET_DOWN_dist2_norm <- TRIPLET_DOWN_dist2/genomewide_ebox_nums$pm_hits[match(names(TRIPLET_DOWN_dist2), genomewide_ebox_nums$ebox)][col(TRIPLET_DOWN_dist2)]
TRIPLET_DOWN_dist2_norm <- TRIPLET_DOWN_dist2_norm %>% tibble::rownames_to_column(var = "bin")
TRIPLET_DOWN_dist2_norm$bin <- as.numeric(TRIPLET_DOWN_dist2_norm$bin)

## plot

all_ebox_ids <- colnames(TRIPLET_DOWN_dist2_norm)[-1]

pdf("../figures/TRIPLET_DEG_EBOX_DISTRIBUTION.pdf", height = 6, width = 8)

for(i in 1:length(all_ebox_ids)){
  
  my_eb <- all_ebox_ids[i]
  
  my_plt <- return_ebox_overlapplot3(up_df = TRIPLET_UP_dist2_norm, down_df = TRIPLET_DOWN_dist2_norm, ebox_id = my_eb) 
  
  print(my_plt)
  
}

dev.off()

## normal theme for MS
pdf("../figures/TRIPLET_DEG_EBOX_DISTRIBUTION_ms.pdf", height = 6, width = 8)

for(i in 1:length(all_ebox_ids)){
  
  my_eb <- all_ebox_ids[i]
  
  my_plt <- return_ebox_overlapplot4(up_df = TRIPLET_UP_dist2_norm, down_df = TRIPLET_DOWN_dist2_norm, ebox_id = my_eb)+
    theme_classic(base_size = 18)+
    theme(axis.ticks.length = unit(0.3, "cm"),
          axis.text = element_text(color = "black"),
          legend.title = element_blank())
  
  print(my_plt)
  
}

dev.off()

## ranking eboxes
## eboxes were ranked according to Myc density (Lin's appraoch)
## then made to heatmap

dmso <- readRDS("../data/HDMB03_DMSO_Ebox_ranking.rds")
ms275 <- readRDS("../data/HDMB03_MS275_Ebox_ranking.rds")

ebox_mat <- dmso %>% 
  dplyr::select(-c(OCCURENCES)) %>% 
  dplyr::rename(DMSO = AVG_HEIGHT) %>% 
  dplyr::left_join(ms275) %>% 
  dplyr::select(-c(OCCURENCES)) %>% 
  dplyr::rename(MS275 = AVG_HEIGHT) %>% 
  dplyr::arrange(desc(DMSO)) %>% 
  dplyr::mutate(diff = DMSO - MS275) %>% 
  tibble::column_to_rownames(var = "EBOX")

ebox_mat_totalEboxNums <- ebox_mat %>% 
  tibble::rownames_to_column(var = "EBOX") %>% 
  dplyr::select(EBOX, DMSO) %>% 
  dplyr::mutate(DMSO = DMSO/max(DMSO)) %>% 
  merge(genomewide_ebox_nums, by.x = "EBOX", by.y = "ebox") %>% 
  dplyr::arrange(desc(DMSO)) %>% 
  tibble::column_to_rownames(var = "EBOX")

ebox_bar_plt <- rowAnnotation(`Genomewide Ebox number (x10^6)` = anno_barplot(ebox_mat_totalEboxNums$pm_hits,
                                                            gp = gpar(fill = "black"),
                                                            axis_param = list(gp = gpar(fontsize = 14),
                                                                              side = "bottom"),
                                                            width = unit(6, "cm")))

ht1_genomewide_ebox <- Heatmap(ebox_mat_totalEboxNums[-2], name = "Myc density",
                               cluster_rows = FALSE,
                               cluster_columns = FALSE,
                               col = colorRampPalette((brewer.pal(n = 7, name = "Reds")))(250),
                               row_names_side = "left",
                               show_column_names = TRUE,
                               row_names_gp = gpar(fontsize = 14, fontface="bold"),
                               column_names_gp = gpar(fontsize = 14, fontface="bold"),
                               column_names_rot = 35,
                               right_annotation = ebox_bar_plt,
                               heatmap_legend_param = list(title = "Myc density", title_gp = gpar(fontsize = 12),
                                                           labels_gp = gpar(fontsize = 14),
                                                           legend_height = unit(4, "cm")))


cell_fun = function(j, i, x, y, width, height, fill) {
  grid.rect(x = x, y = y, width = width * 1, height = height*1, 
            gp = gpar(col = NA, fill = fill))
}


pdf("../figures/DMSO_EboxRanking_onMycDensity_TotalEboxesBar_perMillion.pdf", height = 5, width = 5.5)

ht1_genomewide_ebox
decorate_heatmap_body("Myc density", {
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))
})

dev.off()


pdf("../figures/DMSO_EboxRanking_onMycDensity_TotalEboxesBar_perMillion_2.pdf", height = 5, width = 5)

Heatmap(as.matrix(ebox_mat_totalEboxNums[-2]), name = "Myc density",
        rect_gp = gpar(type = "none"), cell_fun = cell_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        col = colorRampPalette((brewer.pal(n = 7, name = "Reds")))(250),
        row_names_side = "left",
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 14, fontface="bold"),
        column_names_gp = gpar(fontsize = 14, fontface="bold"),
        column_names_rot = 35,
        right_annotation = ebox_bar_plt,
        heatmap_legend_param = list(title = "Myc density", title_gp = gpar(fontsize = 12),
                                    labels_gp = gpar(fontsize = 14),
                                    legend_height = unit(4, "cm")))
decorate_heatmap_body("Myc density", {
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))
})

dev.off()


## public data ebox analysis
## take gene sets where up & down sets are available

all_ebox_up_files <- readRDS("../data/MsigDB_MYC_Targets_up_genesets.rds")

## function to plot up/down ebox dist in one plot
return_up_down_ebox_plot <- function(up_file){
  
  ebox_output <- "../figures/MSigDB_EboxPlots"
  
  down_file = gsub("_UP_10bp_tss", "_DN_10bp_tss", up_file)
  
  ## gene numbers files
  up_num_file <- gsub("_eboxDensities.txt", "", up_file)
  down_num_file <- gsub("_eboxDensities.txt", "", down_file)
  
  up_dat_num <- read.delim(up_num_file, header = FALSE) %>% 
    dplyr::pull(V5) %>% as.character() %>% length()
  
  down_dat_num <- read.delim(down_num_file, header = FALSE) %>% 
    dplyr::pull(V5) %>% as.character() %>% length()
  
  ## title
  my_title <- gsub("_UP_10bp_tss.bed_eboxDensities.txt", "", basename(up_file))
  
  
  ## read file
  readHomer_metagene <- function(file_path){
    
    d = read.delim(file_path, header = TRUE)
    
    cl_names = colnames(d)
    cl_names[1] = "bin"
    
    colnames(d) = cl_names
    
    ## extract bin, *.total.sites* columns
    
    d = d %>% dplyr::select(matches("bin|.total.sites"))
    
    colnames(d) = d %>% colnames() %>% {gsub("\\.total.sites", "", .)}
    
    return(d)
    
  }
  
  ## normalize to the total number of ebox hits genome-wide
  up_raw <- readHomer_metagene(up_file)
  colnames(up_raw) <- gsub("EBOX_", "", colnames(up_raw))
  
  down_raw <- readHomer_metagene(down_file)
  colnames(down_raw) <- gsub("EBOX_", "", colnames(down_raw))
  
  ## normalize to background distribution
  up_raw2 <- up_raw %>% tibble::column_to_rownames(var = "bin")
  down_raw2 <- down_raw %>% tibble::column_to_rownames(var = "bin")
  
  up_raw2_norm <- up_raw2/genomewide_ebox_nums$pm_hits[match(names(up_raw2), genomewide_ebox_nums$ebox)][col(up_raw2)]
  down_raw2_norm <- down_raw2/genomewide_ebox_nums$pm_hits[match(names(down_raw2), genomewide_ebox_nums$ebox)][col(down_raw2)]
  
  up_raw2_norm <- up_raw2_norm %>% tibble::rownames_to_column(var = "bin")
  down_raw2_norm <- down_raw2_norm %>% tibble::rownames_to_column(var = "bin")
  
  ## melt up/down dfs and add extra column explaining group
  ## add Ns for up.down class
  down_lab <- paste0("down_", down_dat_num)
  up_lab <- paste0("up_", up_dat_num)
  
  up_mlt <- up_raw2_norm %>% 
    dplyr::mutate(bin = as.factor(bin)) %>% 
    reshape2::melt() %>% 
    dplyr::mutate(group = "up") 
  
  down_mlt <- down_raw2_norm %>% 
    dplyr::mutate(bin = as.factor(bin)) %>% 
    reshape2::melt() %>% 
    dplyr::mutate(group = "down") 
  
  all_mlt <- rbind(up_mlt, down_mlt)
  all_mlt$bin = as.numeric(levels(all_mlt$bin))[all_mlt$bin]
  
  all_eboxes = all_mlt$variable %>% as.character() %>% unique()
  
  pdf(paste0(ebox_output, "/", my_title, "_UP_DOWN_EBOX_DISTRIBUTION_ms.pdf"), height = 5, width = 7.5)
  
  for(i in 1:length(all_eboxes)){
    
    my_eb <- all_eboxes[i]
    
    my_eb_plt = all_mlt %>% 
      dplyr::filter(variable == get("my_eb")) %>% 
      ggplot(aes(bin, value, color = group))+
      geom_smooth(method = "loess", se = FALSE, size = 1.4)+
      scale_color_manual(values = c("blue", "red"))+
      #theme_vt2()+
      theme_classic(base_size = 18)+
      scale_x_continuous(breaks = c(-2000, -1000, 0, 1000, 2000), labels = c("-2000", "-1000", "TSS", "1000", "2000"))+
      theme(panel.grid = element_blank(),
            legend.title = element_blank(),
            axis.text = element_text(color = "black"),
            axis.ticks.length = unit(0.3, "cm"))+
      ylab("motifs per bp per promoter")+
      xlab("distance from TSS")+
      ggtitle(paste0(my_eb, "; Up-", up_dat_num, ", Down-", down_dat_num))
    
    print(my_eb_plt)
    
  }
  
  dev.off()
  
  
  message(" ## DONE: ", my_title)
  
}

for(i in 1:length(all_ebox_up_files)){
  
  return_up_down_ebox_plot(up_file = all_ebox_up_files[i])
  
}

## Hallmark MYC targtes V2 plots

retunr_ebox_plots <- function(ebox_file_path){
  
  ebox_output <- "../figures/MSigDB_EboxPlots"
  
  ## files and titles
  gene_num_file = gsub("_eboxDensities.txt", "", ebox_file_path)
  gene_num_dat = read.delim(gene_num_file, header = FALSE) %>% 
    dplyr::pull(V5) %>% as.character() %>% length()
  
  gs_title = gsub("_50bp_tss.bed", "", basename(gene_num_file))
  
  ## read homer metagenge files function
  readHomer_metagene <- function(file_path){
    
    d = read.delim(file_path, header = TRUE)
    
    cl_names = colnames(d)
    cl_names[1] = "bin"
    
    colnames(d) = cl_names
    
    ## extract bin, *.total.sites* columns
    
    d = d %>% dplyr::select(matches("bin|.total.sites"))
    
    colnames(d) = d %>% colnames() %>% {gsub("\\.total.sites", "", .)}
    
    return(d)
    
  }
  
  ## wrap fixed scales
  ebox_metageneplot0 <- function(homer_df){
    
    ebox_meta = homer_df
    
    ebox_meta$bin = as.factor(ebox_meta$bin)
    
    ebox_meta_mlt = ebox_meta %>% reshape2::melt()
    
    #as.numeric(levels(f))[f]
    
    ebox_meta_mlt$bin = as.numeric(levels(ebox_meta_mlt$bin))[ebox_meta_mlt$bin]
    
    my_p = ggplot(ebox_meta_mlt, aes(bin, value, color = variable))+
      geom_line(size = 1.5, alpha = 0.4)+
      geom_smooth(method = "loess", se = FALSE, size = 1.5)+
      facet_wrap(. ~ variable, ncol = 4)+
      #theme_bw(base_size = 14)+
      theme_vt_ms()+
      theme(legend.position = "none")+
      ylab("motifs per bp per promoter")+
      xlab("")
    
    return(my_p)
  }
  
  ## wrap fixed scales, only loess
  ebox_metageneplot00 <- function(homer_df){
    
    ebox_meta = homer_df
    
    ebox_meta$bin = as.factor(ebox_meta$bin)
    
    ebox_meta_mlt = ebox_meta %>% reshape2::melt()
    
    #as.numeric(levels(f))[f]
    
    ebox_meta_mlt$bin = as.numeric(levels(ebox_meta_mlt$bin))[ebox_meta_mlt$bin]
    
    my_p = ggplot(ebox_meta_mlt, aes(bin, value, color = variable))+
      #geom_line(size = 1.5, alpha = 0.4)+
      geom_smooth(method = "loess", se = FALSE, size = 1.5)+
      facet_wrap(. ~ variable, ncol = 4)+
      theme_bw(base_size = 14)+
      theme(legend.position = "none")+
      ylab("motifs per bp per promoter")+
      xlab("")
    
    return(my_p)
  }
  
  ## wrap with free scales
  ebox_metageneplot1 <- function(homer_df){
    
    ebox_meta = homer_df
    
    ebox_meta$bin = as.factor(ebox_meta$bin)
    
    ebox_meta_mlt = ebox_meta %>% reshape2::melt()
    
    #as.numeric(levels(f))[f]
    
    ebox_meta_mlt$bin = as.numeric(levels(ebox_meta_mlt$bin))[ebox_meta_mlt$bin]
    
    my_p = ggplot(ebox_meta_mlt, aes(bin, value, color = variable))+
      geom_line(size = 1.5, alpha = 0.4)+
      geom_smooth(method = "loess", se = FALSE, size = 1.5)+
      facet_wrap(. ~ variable, ncol = 4, scales = "free")+
      theme_bw(base_size = 14)+
      theme(legend.position = "none")+
      ylab("motifs per bp per promoter")+
      xlab("")
    
    return(my_p)
  }
  
  ## only loess
  ebox_metageneplot2 <- function(homer_df){
    
    ebox_meta = homer_df
    
    ebox_meta$bin = as.factor(ebox_meta$bin)
    
    ebox_meta_mlt = ebox_meta %>% reshape2::melt()
    
    #as.numeric(levels(f))[f]
    
    ebox_meta_mlt$bin = as.numeric(levels(ebox_meta_mlt$bin))[ebox_meta_mlt$bin]
    
    my_p = ggplot(ebox_meta_mlt, aes(bin, value, color = variable))+
      #geom_line(size = 1.5, alpha = 0.4)+
      geom_smooth(method = "loess", se = FALSE, size = 1.5)+
      facet_wrap(. ~ variable, ncol = 4, scales = "free_y")+
      theme_bw(base_size = 14)+
      theme(legend.position = "none")+
      ylab("motifs per bp per promoter")+
      xlab("")
    
    return(my_p)
  }
  
  ## correcting for genome-wide distribution
  raw_hits <- readHomer_metagene(file_path = ebox_file_path)
  colnames(raw_hits) <- gsub("EBOX_", "", colnames(raw_hits))
  
  ## tip: https://stackoverflow.com/questions/35407852/multiply-columns-with-rows-by-matching-column-name-and-row-name-in-r
  ## data*ref$Values[match(names(data), ref$Names)][col(data)]
  
  ## corrected for background/total occurence of each ebox genomewide
  
  raw_hits2 <- raw_hits %>% tibble::column_to_rownames(var = "bin") 
  
  raw_hits2_norm <- raw_hits2/genomewide_ebox_nums$pm_hits[match(names(raw_hits2), genomewide_ebox_nums$ebox)][col(raw_hits2)]
  raw_hits2_norm <- raw_hits2_norm %>% tibble::rownames_to_column(var = "bin")
  
  ## go with loess, free scales
  
  allEbox_plt = ebox_metageneplot2(raw_hits2_norm)+
    theme_classic(base_size = 18)+
    ggtitle(paste0(gs_title, ", n=", gene_num_dat))+
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.ticks.length = unit(0.3, "cm"),
          panel.grid = element_blank(),
          strip.text.x = element_text(size=15, face = "bold"))
  
  ## only interesting Eboxes from my analysis, CACGTG & CACCTG
  myEbox_plt = raw_hits2_norm %>% 
    dplyr::select(matches("bin|CACGTG|CACCTG")) %>% 
    ebox_metageneplot2()+
    theme_vt()+
    xlab("distance from TSS")+
    scale_color_manual(values = c("black", "black"))+
    ggtitle(paste0(gs_title, ", n=", gene_num_dat))+
    theme(legend.position = "none")
  
  myEbox_CACGTG <- raw_hits2_norm %>% 
    dplyr::select(matches("bin|CACGTG")) %>% 
    dplyr::mutate(bin = as.numeric(bin)) %>% 
    ggplot(aes(x = bin, y = CACGTG))+
    geom_smooth(method = "loess", se = FALSE, size = 1.4, color = "black")+
    theme_classic(base_size = 18)+
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.ticks.length = unit(0.3, "cm"),
          panel.grid = element_blank(),
          strip.text.x = element_text(size=15, face = "bold"))+
    ylab("motifs per bp per promoter")+
    xlab("distance from TSS")+ggtitle("CACGTG")+
    scale_x_continuous(breaks = c(-2000, -1000, 0, 1000, 2000), labels = c("-2000", "-1000", "TSS", "1000", "2000"))
  
  myEbox_CACCTG <- raw_hits2_norm %>% 
    dplyr::select(matches("bin|CACCTG")) %>% 
    dplyr::mutate(bin = as.numeric(bin)) %>% 
    ggplot(aes(x = bin, y = CACCTG))+
    geom_smooth(method = "loess", se = FALSE, size = 1.4, color = "black")+
    theme_classic(base_size = 18)+
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.ticks.length = unit(0.3, "cm"),
          panel.grid = element_blank(),
          strip.text.x = element_text(size=15, face = "bold"))+
    ylab("motifs per bp per promoter")+
    xlab("distance from TSS")+ggtitle("CACCTG")+
    scale_x_continuous(breaks = c(-2000, -1000, 0, 1000, 2000), labels = c("-2000", "-1000", "TSS", "1000", "2000"))
  
  pdf(paste0(ebox_output, "/", gs_title, "_allEboxPlots_ms.pdf"), height = 14, width = 18)
  print(allEbox_plt)
  dev.off()
  
  ## my interesting ebox plt
  pdf(paste0(ebox_output, "/", gs_title, "_CACGTG_CACCTG_ms.pdf"), height = 5, width = 10)
  print(myEbox_plt)
  dev.off()
  
  pdf(paste0(ebox_output, "/", gs_title, "_CACGTG_CACCTG_individual_changedTheme_ms.pdf"), height = 5, width = 6.5)
  print(myEbox_CACGTG)
  print(myEbox_CACCTG)
  dev.off()
  
}

hallmark_myc_v2_tss_ebox_file <- "../data/HALLMARK_MYC_V2_10bp_tss.bed_eboxDensities.txt"

retunr_ebox_plots(ebox_file_path = hallmark_myc_v2_tss_ebox_file)

##


