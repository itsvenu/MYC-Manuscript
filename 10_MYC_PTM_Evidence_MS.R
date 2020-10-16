
setwd("/Users/thatikon/Desktop/ServerView/MYC/scripts/MycHdaci_code/data")

library(tidyverse)

dat <- readxl::read_xlsx("190807_GS_PTMs_evidence.xlsx")

pdf("../figures/MYC_PTM_Evidence_table.pdf", height = 10, width = 35)
ggpubr::ggtexttable(dat, rows = NULL)
dev.off()

pdf("../figures/MYC_PTM_Evidence_table_01.pdf", height = 10, width = 20)

dat %>% 
  dplyr::select(c(Position:X__6)) %>% 
  ggpubr::ggtexttable(rows = NULL)

dev.off()

##
pdf("../figures/MYC_PTM_Evidence_table_02.pdf", height = 10, width = 20)
dat %>% 
  dplyr::select(Position:`Modified sequence`, `Log10 Intensity`:X__11) %>% 
  ggpubr::ggtexttable(rows = NULL)

dev.off()