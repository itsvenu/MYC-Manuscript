######################################
## Venu Thatikonda                  ##
## v.thatikonda@kitz-heidelberg.de  ##
## @nerd_yie                        ##
######################################

## source: /b06x-lsdf/infra5-lsdf/mbHDAC2/scripts/DRC_Rescue_Experiment/DRC_Rescue_plots.R

setwd("~/Desktop/ServerView/MYC/scripts/MycHdaci_code/scripts")

library(tidyverse)
library(cowplot)

## custom theme
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

theme_vt_cl <- function(){
  theme_classic(base_size=18) %+replace%
    theme(axis.ticks.length = unit(0.3, "cm"),
          #panel.border = element_rect(size=1.3, color="black", fill=NA),
          axis.text = element_text(color="black"),
          axis.title = element_text(color="black"),
          legend.text = element_text(color="black"),
          legend.title = element_text(color="black")
    )
}

###################################
###### 0 Rescue experiment ########
###################################

# data
# INPUT <- "~/Desktop/ServerView/MYC/results/DRC_Rescue_Experiment"
# res_rescue <- read.delim(paste0(INPUT, "/Results_Rescue_vt.txt"), header = FALSE)
# saveRDS(res_rescue, file = "data/Results_Rescue_vt.rds")

res_rescue <- readRDS("data/Results_Rescue_vt.rds")

# no Dox/Myc low <- noDox_MycLow
# Dox/Myc High <- Dox_MycHigh
# Dox + siMyc/low Myc <- Dox_siMyc_MycLow

colnames(res_rescue) <- c("condition",
                          "noDox_MycLow_R1", "noDox_MycLow_R2", "noDox_MycLow_R3", "noDox_MycLow_R4",
                          "Dox_MycHigh_R1", "Dox_MycHigh_R2", "Dox_MycHigh_R3", "Dox_MycHigh_R4",
                          "Dox_siMyc_MycLow_R1", "Dox_siMyc_MycLow_R2", "Dox_siMyc_MycLow_R3", "Dox_siMyc_MycLow_R3")

## calculate Mean-sd and plot barplots
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

## add error bars
res_rescue_eb <- res_rescue %>% reshape2::melt() %>%
  dplyr::mutate(variable = gsub("_R.*", "", variable)) %>% 
  summarySE(measurevar = "value", groupvars = c("condition", "variable"))

## order of features
res_rescue_eb$condition <- factor(res_rescue_eb$condition, levels = c("untreated", "DMSO", "0.2µM", "1µM", "5µM"))

res_rescue_eb$variable <- factor(res_rescue_eb$variable, levels = c("noDox_MycLow", "Dox_MycHigh", "Dox_siMyc_MycLow"))

res_rescue_plt <- ggplot(res_rescue_eb, aes(x = condition, y = value, fill = variable))+
  geom_bar(position = position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin = value-se, ymax = value+se),
                width=.3, size = 1.1,
                position=position_dodge(.9))+
  scale_fill_manual(labels = c("No Dox/Myc low", "Dox/Myc high", "Dox siMyc/Myc low"),
                    values = c("darkgray", "red", "lightblue"))+
  theme_vt()+
  xlab("")+
  ylab("number of viable cells (% of untreated)")+
  scale_y_continuous(labels = scales::percent_format(scale=1))+
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "top")

res_rescue_plt2 <- ggplot(res_rescue_eb, aes(x = condition, y = value, fill = variable))+
  geom_bar(position = position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin = value-se, ymax = value+se),
                width=.3, size = 1.1,
                position=position_dodge(.9))+
  scale_fill_manual(labels = c("No Dox/Myc low", "Dox/Myc high", "Dox siMyc/Myc low"),
                    values = c("darkgray", "red", "lightblue"))+
  theme_vt_cl()+
  xlab("")+
  ylab("number of viable cells (% of untreated)")+
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), 
                     labels = c("0%", "25%", "50%", "75%", "100%"),
                     expand = expand_scale(mult = c(0, 0.1)))+
  scale_x_discrete(expand = expand_scale(add = 0.5))+
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = c(0.8, 0.9),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.line.x.bottom = element_blank(),
        axis.ticks.x = element_blank())

### change figure legends
res_rescue_eb2 <- res_rescue_eb %>% 
  dplyr::mutate(new_legend = dplyr::case_when(variable == "noDox_MycLow" ~ "no dox - MYC low",
                                              variable == "Dox_MycHigh" ~ "+dox - MYC high",
                                              variable == "Dox_siMyc_MycLow" ~ "+dox + siMYC - MYC low")) 


res_rescue_plt3 <- ggplot(res_rescue_eb2, aes(x = condition, y = value, fill = new_legend))+
  geom_bar(position = position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin = value-se, ymax = value+se),
                width=.3, size = 1.1,
                position=position_dodge(.9))+
  scale_fill_manual(labels = c("no dox - MYC low", "+dox - MYC high", "+dox + siMYC - MYC low"),
                    values = c("darkgray", "red", "lightblue"))+
  theme_vt_cl()+
  xlab("")+
  ylab("number of viable cells (% of untreated)")+
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), 
                     labels = c("0%", "25%", "50%", "75%", "100%"),
                     expand = expand_scale(mult = c(0, 0.1)))+
  scale_x_discrete(expand = expand_scale(add = 0.5))+
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = c(0.8, 0.9),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.ticks.x = element_blank())


save_plot("../figures/Results_rescue_plot_vt.pdf", res_rescue_plt2, base_aspect_ratio = 1.3, base_height = 6, base_width = 8)

save_plot("../figures/Results_rescue_plot_vt2.pdf", res_rescue_plt3, base_aspect_ratio = 1.3, base_height = 6, base_width = 8)

# pdf(paste0(INPUT, "/Results_rescue_plot_vt.pdf"), height = 7, width = 10)
# print(res_rescue_plt)
# dev.off()

pdf(paste0("../figures", "/Results_rescue_plot_vt.pdf"), height = 6, width = 8)
print(res_rescue_plt)
dev.off()

