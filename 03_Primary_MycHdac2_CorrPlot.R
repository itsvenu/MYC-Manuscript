## source: /b06x-lsdf/infra5-lsdf/mbHDAC2/scripts/data_290917/Figures/RNApol_Promoters/MYC_HDAC2_lineplot.R

setwd("~/Desktop/ServerView/MYC/scripts/MycHdaci_code/scripts")

library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)

# theme_vt <- function(){
#   theme_classic(base_size=18) %+replace%
#     theme(axis.ticks.length = unit(0.3, "cm"),
#           panel.border = element_rect(size=1, color="black", fill=NA),
#           axis.text = element_text(color="black", face = "bold"),
#           axis.title = element_text(color="black", face = "bold"),
#           legend.text = element_text(color="black", face = "bold"),
#           legend.title = element_text(color="black", face = "bold")
#     )
# }

theme_vt <- function(){
  theme_classic(base_size=18) %+replace%
    theme(axis.ticks.length = unit(0.3, "cm"),
          panel.border = element_rect(size=1.3, color="black", fill=NA),
          axis.text = element_text(color="black"),
          axis.title = element_text(color="black"),
          legend.text = element_text(color="black"),
          legend.title = element_text(color="black")
    )
}

d = readRDS("../data/AllFactors_Enrichments_ProteinCodingPromoters.rds")

d.mean = d %>% dplyr::mutate(PROMOTER = paste(chr, start, end, sep = "_")) %>%
  dplyr::mutate(MYC = (MB1_cMYC + MB2_cMYC + MB3_cMYC)/3,
                HDAC2 = (MB1_HDAC2 + MB2_HDAC2 + MB3_HDAC2)/3,
                RNApol = (MB1_RNApolII + MB2_RNApolII + MB3_RNApolII)/3,
                H3K27ac = (MB1_H3K27ac + MB2_H3K27ac + MB3_H3K27ac)/3) %>%
  dplyr::select(PROMOTER, MYC, HDAC2, RNApol, H3K27ac) %>%
  dplyr::arrange(MYC)

d.mean$PROMOTER = factor(d.mean$PROMOTER, levels = d.mean$PROMOTER)


d.mean.mlt = d.mean %>% dplyr::select(PROMOTER, MYC, HDAC2) %>% reshape2::melt()

# ggline(d.mean.mlt, x = "PROMOTER", y = "value", color = "variable", plot_type = "l")+
#   theme(axis.text.x = element_blank())

## MYC

p0 = ggplot(d.mean, aes(MYC, HDAC2)) +
  geom_point(size = 0.1, alpha=0) +
  geom_smooth() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  xlab("input normalized MYC signal (log2)")+
  ylab("input normalized HDAC2 signal (log2)")+
  theme_vt()

## scatter plot with corelation value

p1 = ggscatter(d.mean, x = "MYC", y = "HDAC2", add = "reg.line", add.params = list(color = "blue"),
               color = "black", size = 0.3, alpha = 0.7) +
  stat_cor(method = "spearman", size = 5) +
  xlab("input normalized MYC signal (Log2)")+
  ylab("input normalized HDAC2 signal (Log2)")+
  theme_vt()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

# 
# p2 <- ggscatter(d.mean, x = "MYC", y = "HDAC2", add = "reg.line", add.params = list(color = "blue"),
#           color = "black", size = 0.3, alpha = 0.7) +
#   stat_cor(method = "spearman", size = 6) +
#   xlab("input normalized MYC signal (Log2)")+
#   ylab("input normalized HDAC2 signal (Log2)")+
#   theme_minimal(base_size = 18)+
#   theme(panel.grid.minor = element_blank(),
#         axis.text = element_text(color = "black"))

p2 <- ggscatter(d.mean, x = "MYC", y = "HDAC2", add = "reg.line", add.params = list(color = "blue"),
          color = "black", size = 0.3, alpha = 0.7) +
  stat_cor(method = "spearman", size = 6) +
  xlab("input normalized MYC signal (Log2)")+
  ylab("input normalized HDAC2 signal (Log2)")+
  theme_classic(base_size = 18)+
  theme(axis.text = element_text(color = "black"),
        axis.ticks.length = unit(0.3, "cm"))


p3 <- ggscatter(d.mean, x = "MYC", y = "HDAC2", add = "reg.line", add.params = list(color = "blue"),
                color = "black", size = 0.3, alpha = 0.7) +
  stat_cor(method = "pearson", size = 6) +
  xlab("input normalized MYC signal (Log2)")+
  ylab("input normalized HDAC2 signal (Log2)")+
  theme_classic(base_size = 18)+
  theme(axis.text = element_text(color = "black"),
        axis.ticks.length = unit(0.3, "cm"))

pdf("../figures/MYC_HDAC2_corr_proteinCodingPromoters.pdf", height = 6, width = 8)
print(p1)
dev.off()

save_plot("../figures/MYC_HDAC2_corr_proteinCodingPromoters.pdf", p1, base_aspect_ratio = 1.2, base_height = 6)

save_plot("../figures/MYC_HDAC2_corr_proteinCodingPromoters_minsize.pdf", p1, base_aspect_ratio = 1.3, base_height = 4, base_width = 4)

save_plot("../figures/MYC_HDAC2_corr_proteinCodingPromoters_minsize_minimal.pdf", p2, base_aspect_ratio = 1.3, base_height = 6)


# 1kb promoters to be consistent across all figures -----------------------

dat_1kb <- readRDS("../data/AllFactors_Enrichments_ProteinCodingPromoters_1kb.rds")

dat_1kb_fmt <- dat_1kb %>% dplyr::mutate(PROMOTER = paste(chr, start, end, sep = "_")) %>%
  dplyr::mutate(MYC = (MB1_cMYC + MB2_cMYC + MB3_cMYC)/3,
                HDAC2 = (MB1_HDAC2 + MB2_HDAC2 + MB3_HDAC2)/3,
                RNApol = (MB1_RNApolII + MB2_RNApolII + MB3_RNApolII)/3,
                H3K27ac = (MB1_H3K27ac + MB2_H3K27ac + MB3_H3K27ac)/3) %>%
  dplyr::select(PROMOTER, MYC, HDAC2, RNApol, H3K27ac) %>%
  dplyr::arrange(MYC) 

p4_1kb <- ggscatter(dat_1kb_fmt, x = "MYC", y = "HDAC2", add = "reg.line", add.params = list(color = "blue"),
          color = "black", size = 0.3, alpha = 0.7) +
  stat_cor(method = "pearson", size = 6) +
  xlab("input normalized MYC signal (Log2)")+
  ylab("input normalized HDAC2 signal (Log2)")+
  theme_classic(base_size = 18)+
  theme(axis.text = element_text(color = "black"),
        axis.ticks.length = unit(0.3, "cm"))

save_plot("../figures/MYC_HDAC2_corr_proteinCodingPromoters1kb_minsize_minimal.pdf", p4_1kb, base_aspect_ratio = 1.3, base_height = 6)

## HALLMARK MYC TARGETS binding signals

promoters_1kb <- read.delim("~/Desktop/ServerView/MYC/results/data_290917/gencode-transcripts/gencode.v19.proteincoding.genelevel.PROMOTER_-1kb_+1kb.bed", header = FALSE)
promoters_1kb <- promoters_1kb %>% 
  dplyr::mutate(PROMOTER = paste(V1, V2, V3, sep = "_")) %>% 
  dplyr::select(PROMOTER, V4) 

MYC_TARGETS <- "~/Desktop/ServerView/MYC/results/data_290917/public-data"

v1 <- read.delim(paste0(MYC_TARGETS, "/MYC_TARGETS_V1.txt"), header = FALSE)
v2 <- read.delim(paste0(MYC_TARGETS, "/MYC_TARGETS_V2.txt"), header = FALSE)

v1v2_common <- intersect(v1$V1, v2$V1) %>% 
  as.data.frame() %>% 
  dplyr::rename(V1 = ".")

## remove these from V1
v1 <- setdiff(v1$V1, v1v2_common$V1) %>% 
  as.data.frame() %>% 
  dplyr::rename(V1 = ".")


g19_signals <- promoters_1kb %>% 
  dplyr::left_join(dat_1kb_fmt) %>% 
  dplyr::filter(MYC != "NA") %>% 
  dplyr::select(-c(PROMOTER))

## version 2 signals
v2_status <- g19_signals %>% 
  merge(v2, by.x = "V4", by.y = "V1") %>% 
  dplyr::select(-c(RNApol)) %>% 
  dplyr::mutate(class = "V2") %>% 
  dplyr::mutate(group = dplyr::case_when(MYC > 0.5 & HDAC2 > 0.5 & H3K27ac > 0.5 ~ "Triplet bound",
                                         MYC <= 0.5 | HDAC2 <= 0.5 | H3K27ac <= 0.5 ~ "Non triplet"))

v2_status_fractions <- v2_status %>% 
  dplyr::count(group) 

v2_status_fractions$fractions <- v2_status_fractions$n/sum(v2_status_fractions$n)
v2_status_fractions$class = "V2"

##
v1_status <- g19_signals %>% 
  merge(v1, by.x = "V4", by.y = "V1") %>% 
  dplyr::select(-c(RNApol)) %>% 
  dplyr::mutate(class = "V1") %>% 
  dplyr::mutate(group = dplyr::case_when(MYC > 0.5 & HDAC2 > 0.5 & H3K27ac > 0.5 ~ "Triplet bound",
                                         MYC <= 0.5 | HDAC2 <= 0.5 | H3K27ac <= 0.5 ~ "Non triplet",
                                         MYC == "NA" | HDAC2 == "NA" | H3K27ac == "NA" ~ "NO-SIGNAL"))


v1_status_fractions <- v1_status %>% 
  dplyr::count(group) 

v1_status_fractions$fractions <- v1_status_fractions$n/sum(v1_status_fractions$n)
v1_status_fractions$class = "V1"

###

v1v2_fractionTriplet_plt <- rbind(v1_status_fractions, v2_status_fractions) %>% 
  ggbarplot(x = "class", y = "fractions", fill = "group", color = NA) +
  scale_fill_manual(values = c("#D3D3D3", "#696969"))+
  theme_vt_cl()+
  xlab("")+
  ylab("% of HALLMARK_MYC_TARGET\ngenes")+
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = c("0%", "25%", "50%", "75%", "100%"),
                     expand = expand_scale(mult = c(0, 0.1)))+
  scale_x_discrete(expand = expand_scale(add = 0.5))+
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.ticks = element_line(color = "black"),
        axis.ticks.x = element_blank(),
        axis.line.x.bottom = element_blank())

pdf("../figures/Primary_HALLMARK_MYC_TARGETS_V1V2_FractionTripletBound.pdf", height = 5, width = 6.5)
print(v1v2_fractionTriplet_plt)
dev.off()






