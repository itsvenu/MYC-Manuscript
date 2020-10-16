
## prepare beta values for cpg_* for 3 primary tumors
## for GEO submission

setwd("~/Desktop/ServerView/MYC/scripts/MycHdaci_code/scripts")

library(dplyr)
library(minfi)

input <- "~/Desktop/ServerView/MYC/scripts/450k_meth/idat"

all_files <- list.files(path = input, 
                        pattern = "*.idat$", full.names = TRUE)

all_files <- gsub("_Grn.idat", "", all_files) %>% 
  gsub("_Red.idat", "", .) %>% unique()

meth_450k <- read.metharray(basenames = all_files, verbose = TRUE, force = TRUE) 

meth_450k_rg <- preprocessRaw(meth_450k)

RSet <- ratioConvert(meth_450k_rg, what = "both", keepCN = TRUE)

betas <- getBeta(RSet)

betas <- betas %>% as.data.frame() %>% 
  dplyr::rename(MB1 = `201490030134_R02C01`,
                MB2 = `201490030134_R01C01`,
                MB3 = `200876170031_R06C01`) %>% 
  tibble::rownames_to_column(var = "probe") 

write.table(betas, file = paste0(input, "/Primary_Betas_450k.txt"), sep = "\t", row.names = F, quote = F)




