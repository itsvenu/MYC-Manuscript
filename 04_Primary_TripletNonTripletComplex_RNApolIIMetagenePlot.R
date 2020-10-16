## source: /b06x-lsdf/infra5-lsdf/mbHDAC2/scripts/affy-ge/all_timepoints_deg/closestEnhancers/metaplot_smoothing.R
## Triplet_complex vs non_complex RNApolII metagene plot

setwd("~/Desktop/ServerView/MYC/scripts/MycHdaci_code/scripts")

library(tidyverse)

# theme_vt2 <- function(){
#   theme_classic(base_size=18) %+replace%
#     theme(axis.ticks.length = unit(0.3, "cm"),
#           panel.border = element_rect(size=1, color="black", fill=NA),
#           axis.text = element_text(color="black", face = "bold"),
#           axis.title = element_text(color="black", face = "bold"),
#           legend.text = element_text(color="black", face = "bold"),
#           legend.title = element_text(color="black", face = "bold")
#     )
# }

theme_vt2 <- function(){
  theme_classic(base_size=18) %+replace%
    theme(axis.ticks.length = unit(0.3, "cm"),
          panel.border = element_rect(size=1.3, color="black", fill=NA),
          axis.text = element_text(color="black"),
          axis.title = element_text(color="black"),
          legend.text = element_text(color="black"),
          legend.title = element_text(color="black")
    )
}

## function to format `deeptools plotProfile` output
return_plotProfile_formated_RNApolII_primaryTriplet <- function(file_path){
  
  # file_path = tripletComplex_nonComplex
  dat <- read.delim(file_path, header = TRUE, skip = 1)
  
  dat_mlt <- dat %>%
    dplyr::select(-c(bins)) %>%
    tibble::column_to_rownames(var = "X") %>% t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "bin") %>%
    dplyr::filter(bin != "X") %>%
    dplyr::mutate(bin = gsub("X", "", bin)) %>%
    dplyr::rename(Triplet = Triplet_complex, Non_triplet = non_complex) %>%
    reshape2::melt() %>%
    dplyr::mutate(bin = as.numeric(bin))
  
  return(dat_mlt)
  
}

## function to plot metagene plot the data returned by above function
plotProfile_plot_RNApolII_themeClassic <- function(melted_dat){
  
  # melted_dat <- up_rnapol_mlt
  
  ## custom labels
  lab_pos <- c(100, 400)
  lab_nm <- c("TSS", "TES")
  
  spline_int <- as.data.frame(spline(melted_dat$bin, melted_dat$value))
  
  my_plt <- ggplot(melted_dat, aes(bin, value, color = variable))+
    geom_line(size = 2)+
    scale_x_continuous(breaks = lab_pos, labels = lab_nm, expand = expand_scale(mult = c(0, 0)))+
    theme_vt2()+
    xlab("")+
    ylab("input normalized RNApolII signal (Log2)")+
    #scale_color_manual(values = c("RNApolII_untreated" = "darkgray", "RNApolII_treated" = "black"))+
    theme(legend.title = element_blank(),
          legend.position = c(0.85, 0.9))
  return(my_plt)
  
}

## analysis
tripletComplex_nonComplex <- "../data/TRIPLETCOMPLEX_NONCOMPLEX_EXPRESSED_RNApolII_grays.data"

tripletComplex_nonComplex_mlt <- return_plotProfile_formated_RNApolII_primaryTriplet(file_path = tripletComplex_nonComplex)

tripletComplex_nonComplex_mlt2 <- tripletComplex_nonComplex_mlt %>% 
  dplyr::mutate(variable = gsub("_", " ", variable))
tripletComplex_nonComplex_mlt2$variable <- factor(tripletComplex_nonComplex_mlt2$variable, levels = c('Triplet', 'Non triplet'))

##
pdf("../figures/TRIPLET_vs_NONTRIPLET_RNApolII_metagene.pdf", height = 5.5, width = 8)

plotProfile_plot_RNApolII_themeClassic(melted_dat = tripletComplex_nonComplex_mlt2)+
  scale_color_manual(values = c("Triplet" = "black", "Non triplet" = "darkgray"))

dev.off()

triplet_nontriplet_metagene_primary <- plotProfile_plot_RNApolII_themeClassic(melted_dat = tripletComplex_nonComplex_mlt2)+
  scale_color_manual(values = c("Triplet" = "black", "Non triplet" = "darkgray"))


triplet_nontriplet_metagene_primary_classic <- plotProfile_plot_RNApolII_themeClassic(melted_dat = tripletComplex_nonComplex_mlt2)+
  scale_color_manual(values = c("Triplet" = "black", "Non triplet" = "darkgray"))+
  theme_classic(base_size = 18)+
  scale_y_continuous(breaks = c(0, 1, 2, 3), labels = c("0", "1", "2", "3"))+
  theme(axis.text = element_text(color = "black"),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_blank(),
        legend.position = c(0.7, 0.9))

save_plot("../figures/TRIPLET_vs_NONTRIPLET_RNApolII_metagene.pdf",triplet_nontriplet_metagene_primary, base_aspect_ratio = 1.3, base_height = 5.5, base_width = 8)

save_plot("../figures/TRIPLET_vs_NONTRIPLET_RNApolII_metagene_classic.pdf", triplet_nontriplet_metagene_primary_classic, base_aspect_ratio = 1.3, base_height = 5, base_width = 7.5)

