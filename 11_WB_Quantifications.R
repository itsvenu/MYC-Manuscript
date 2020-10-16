## WB quantiufications

library(tidyverse)

dat <- readxl::read_xlsx("../data/wb_quantification.xlsx")

dat <- dat %>% 
  reshape2::melt() 

dat$rep <- c(rep("Un", 3), rep("Tr", 3),
             rep("Un", 3), rep("Tr", 3))

dat$gg <- c(rep(c("HA1", "HA2", "HA3"), 2), rep(c("Kless1", "Kless2", "Kless3"), 2))

dat$rep <- factor(dat$rep, levels = c("Un", "Tr"))

# df <- with(df, df[order(gram, -as.numeric(number)), ])

# dat2 <- with(dat, dat[order(rep, -as.numeric(rep)), ])

dat$variable <- factor(dat$variable, levels = c("HA-untreated", "HA-treated",
                                                "Kless-untreated", "Kless-treated"))

# ggplot(dat, aes(x = gg, y = value, fill = rep))+
#    geom_bar(stat = "identity", position="dodge")

wb_plt <- ggbarplot(dat, x = "gg", y = "value", fill = "rep", color = NA,
          position = position_dodge2(preserve = "single"))+
  scale_fill_manual(values = c("gray", "red"))+
  theme_vt_cl()+
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)))+
  scale_x_discrete(expand = expand_scale(add = 0.5))+
  ylab("MYC intensity to HDAC2 control")+
  xlab("")+
  theme(axis.ticks = element_line(color = "black"),
        axis.ticks.x = element_blank(),
        axis.line.x.bottom = element_blank())

pdf("../figures/WB_HA_Kless_quantifications.pdf", height = 5, width = 7)
print(wb_plt)
dev.off()

## mean and error bars plot
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

dat_ci <- dat %>% 
  dplyr::select(-c(gg)) %>% 
  summarySE(measurevar = "value", groupvars = c("variable", "rep"))

dat_ci$gg <- c("HA", "HA", "Kless", "Kless")

wb_ci_plt <- ggplot(dat_ci, aes(x = gg, y = value, fill = rep))+
  geom_bar(position = position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin = value-se, ymax = value+se),
                width=.3, size = 1.1,
                position=position_dodge(.9))+
  scale_fill_manual(values = c("gray", "red"))+
  theme_vt_cl()+
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)),
                     breaks = c(0, 1, 2, 2.5),
                     labels = c(0, 1, 2, 2.5))+
  scale_x_discrete(expand = expand_scale(add = 0.5))+
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = c(0.8, 0.9),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.line.x.bottom = element_blank(),
        axis.ticks.x = element_blank())+
  ylab("MYC intensity to HDAC2 control")+
  xlab("")
  
pdf("../figures/WB_HA_Kless_quantifications_errorBars.pdf", height = 5, width = 4)
print(wb_ci_plt)
dev.off()

## Wb quantification on chromatin
wb_c <- readxl::read_xlsx("../data/wb_quantification_chromatin.xlsx")

wb_c_confI <- wb_c %>% 
  reshape2::melt() %>% 
  summarySE(measurevar = "value", groupvars = c("variable"))

wb_c_confI_plt <- ggplot(wb_c_confI, aes(x = variable, y = value, fill = variable))+
  geom_bar(position = position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin = value-se, ymax = value+se),
                width=.3, size = 1.1,
                position=position_dodge(.9))+
  scale_fill_manual(values = c("gray", "red"))+
  theme_vt_cl()+
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)))+
  scale_x_discrete(expand = expand_scale(add = 0.5))+
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.line.x.bottom = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  ylab("MYC intensity to HDAC2 control")+
  xlab("")

pdf("../figures/WB_MYC_on_chromatin_quantifications.pdf", height = 4, width = 2.5)
print(wb_c_confI_plt)
dev.off()


