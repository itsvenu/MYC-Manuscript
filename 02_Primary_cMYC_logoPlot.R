## source: /b06x-lsdf/infra5-lsdf/mbHDAC2/scripts/data_290917/Figures/cMYC_logo/cMYC_logo.R

## cMYC logo-plot from homer

setwd("~/Desktop/ServerView/MYC/scripts/MycHdaci_code/scripts")

library(tidyverse)
library(ggplot2)
library(ggseqlogo)
library(ggpubr)
library(gridExtra)
library(grid)

# theme_vt <- function(){
#   theme_classic(base_size=18) %+replace%
#     theme(axis.ticks.length = unit(0.3, "cm"),
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

## cmyc pfm
cmyc <- readRDS("../data/cMYC_PFM_homer.rds")

cmyc_mat <- cmyc %>% t()
rownames(cmyc_mat) <- c("A", "C", "G", "T")

logo_plt <- ggseqlogo(cmyc_mat)+
  theme_vt()+
  xlab("Position")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())

## gtable pvalues
## cMYC pvalues from homer motif analysis

df <- tibble(`MYC Peaks` = "1e-125",
                 `HDAC2 Peaks` = "1e-53")

rownames(df) <- c("P.value")

pval_table <- ggtexttable(df %>% as.data.frame(), theme = ttheme("mGreen", base_size = 34))

pdf("../figures/cMYCLogo_Pvals.pdf", height = 3, width = 7.7)
print(logo_plt)
print(pval_table)
dev.off()

save_plot("../figures/cMYCLogo_Pvals.pdf", logo_plt, base_aspect_ratio = 1.3, base_height = 3, base_width = 5)

### motif plot from dreme results

# 0.000000 1.000000 0.000000 0.000000
# 0.000000 1.000000 0.000000 0.000000
# 0.000000 1.000000 0.000000 0.000000
# 0.000000 1.000000 0.000000 0.000000
# 0.169960 0.000000 0.324111 0.505929
# 0.000000 1.000000 0.000000 0.000000
# 0.000000 1.000000 0.000000 0.000000
# 0.000000 1.000000 0.000000 0.000000

dreme_motif <- read.delim("../data/CCCCDCCC_freqs.txt", header = FALSE, sep = " ")
colnames(dreme_motif) <- c("A", "C", "G", "T")

dreme_motif_mat <- dreme_motif %>% t()

sp1_plt <- ggseqlogo(dreme_motif_mat)+
  xlab("Position")+
  theme_classic(base_size = 18)+
  theme(axis.text = element_text(color = "black"),
        axis.ticks.length = unit(0.3, "cm"))+
  ggtitle("E-value: 2.7e-004")

save_plot("../figures/_supp/DREME_SP1_motifLogo.pdf", sp1_plt, base_aspect_ratio = 1.3, base_height = 3, base_width = 5)

### de novo result plots ####

## from HDAC2 peaks, best e-box hit
hdac2_myc_denovo <- "/Users/thatikon/Desktop/ServerView/MYC/results/data_290917/Figures/HDAC2_motif_on_MYC/motifs_top3k_summits/HDAC2_3K_denovo/homerResults/motif2.similar1.motif"
myc_myc_denovo <- "~/Desktop/ServerView/MYC/results/data_290917/Figures/HDAC2_motif_on_MYC/motifs_top3k_summits/MYC_3K_denovo/homerResults/motif3.motif"

process_homer <- function(x){
  
  x_rd <- data.table::fread(x, skip = 1) %>%
    dplyr::rename(A = "V1", C = "V2", G = "V3", T = "V4")
  
  return(x_rd)
}

hdac2_myc_denovo <- process_homer(hdac2_myc_denovo) %>% t()
myc_myc_denovo <- process_homer(myc_myc_denovo) %>% t()

pdf("../figures/HDAC2_denovo_MYCmotif.pdf", height = 3, width = 5)

ggseqlogo(hdac2_myc_denovo)+
  theme_classic(base_size = 18)+
  xlab("Position")+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.ticks.length = unit(0.3, "cm"))+
  ggtitle("Pvalue = 1e-82")

dev.off()

pdf("../figures/MYC_denovo_MYCmotif.pdf", height = 3, width = 5)

ggseqlogo(myc_myc_denovo)+
  theme_classic(base_size = 18)+
  xlab("Position")+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.ticks.length = unit(0.3, "cm"))+
  ggtitle("Pvalue = 1e-124")

dev.off()






