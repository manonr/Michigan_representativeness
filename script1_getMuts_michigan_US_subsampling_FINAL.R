# ============================================================================
# script1_getMuts_michigan_US_subsampling_FINAL.R
# Step 1: Identify high-frequency mutation sites and define haplotypes.
# For each subtype/threshold/season combination, identifies variable AA
# positions, constructs haplotype frequency tables (with/without Michigan),
# and computes Shannon entropy.
#
# Inputs: H1_ALLmuts.rds, H3_ALLmuts.rds
# Outputs: haploDF_*.csv, bigResultsDF_*.csv
# Dependencies: functions_pick_strains_functions.R
# ============================================================================


library(tidyr)
library(dplyr)
require(ggplot2)
require(ggseqlogo)
library(bioseq)
library(ape)
library(stringr)
library(ggpubr)
library(readr)
library(lubridate)
library(reshape2)
library(lubridate)
library(tibble)
library(readr)
library(Biostrings)

# First, use auspice to convert json to dataframe
# # # in terminal
#  nextstrain shell .
 # python3 auspice_tree_to_table.py \
 # --tree auspice/HA_H1_nextstrain.json \
 # --output-metadata HA_H1_nextstrain.json.tsv \
 # --output-tree HA_H1_nextstrain.tsv.nwk \
 # --include-internal-nodes

rootDir <- paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/")
setwd(rootDir)

source("/Users/manonragonnet/OneDrive/Github/Michigan_rep/R_Michigan/Representativeness/functions_pick_strains_functions.R")

subtypes <- c("H1", "H3")

#for (subtype in subtypes){            ###should not need to be rerun
#  USTable <- create_big_mut_df(subtype, "runs3")
#  saveRDS(USTable, file=paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/", subtype,"_ALLmuts.rds"))
#}


thresholds <- c( 0.2, 0.15,0.1, 0.05, 0.01, 0.005)#0.05, 0.1, 0.2)
mySeasons <- c(2021, 2022, 2023, 2024, paste(2021:2024, collapse="-"))

# first create the big dataframes
# then straight version, no reps/ drops

resultsDF <- data.frame(subtype = "H1", threshold = 0.05, mySeasons = 2021,
                        nSeq = 0, nMichSeq= 0,nVarSites= 0, nhighFreqSites= 0, 
                        nHaplo_wMich = 0, nHaplo_noMich = 0, highfreqMutSites="none",
                        UMichEntropy=0, nonMichEntropy = 0, allEntropy=0)
i <- 0
for (subtype in subtypes){
  
  #read once
  bigMutDF <- readRDS(paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/", subtype,"_ALLmuts.rds"))
  cat("[script1] ", subtype, " | loaded RDS | total:", nrow(bigMutDF),
      "| UMich:", sum(bigMutDF$UMich == "UMich"), "| notUOM:", sum(bigMutDF$UMich == "notUOM"), "\n")

  for (threshold in thresholds){
    for (mySeason in mySeasons){
      print(paste(subtype, threshold, mySeason))
      if(length(strsplit(as.character(mySeason), "-")[[1]])>1){
        mySeason <- base::strsplit(as.character(mySeason), "-")[[1]]
      }
      
      # read RDS
      seasonMutDF <- bigMutDF %>% filter(season %in% mySeason)
      allMutatingSites <- length((getMutList(seasonMutDF, "aapos", myCutoff=1/nrow(seasonMutDF))))
     
      runName <- paste(subtype, threshold, paste(mySeason, collapse = "."), sep="_")
      
      seasonMutDF2 <- create_mut_df(seasonMutDF,threshold = threshold, mySeason = mySeason)
     
      
      if(!is.null(seasonMutDF2)){
        cat("[script1] ", subtype, " | threshold", threshold, "| season", paste(mySeason, collapse="."),
            "| total:", nrow(seasonMutDF2),
            "| UMich:", sum(seasonMutDF2$UMich == "UMich"),
            "| notUOM:", sum(seasonMutDF2$UMich == "notUOM"), "\n")
        i<- i+1
        haploDF_noMich <- haploOrder(seasonMutDF2, inclMich = F)
        haploDF_wMich <- haploOrder(seasonMutDF2, inclMich = T)
        haploDF_MichOnly <- haploOrder(seasonMutDF2, inclMich = "only")
        write.csv(haploDF_noMich, paste0(rootDir, "haplo_stats3/haploDF_noMich_", runName,".csv"), row.names = F)
        write.csv(haploDF_wMich, paste0(rootDir, "haplo_stats3/haploDF_wMich_", runName,".csv"), row.names = F)
        
        #entropy
        UMichEntropy <- round(shannon_entropy(haploDF_MichOnly$Count),3)
        nonMichEntropy <- round(shannon_entropy(haploDF_noMich$Count),3)
        allEntropy <- round(shannon_entropy(haploDF_wMich$Count),3)
        print(allEntropy)
        
        # stats and results table
        resultsDF[i,] <- c(subtype, threshold, paste(mySeason, collapse = "."),
                           nrow(seasonMutDF2), nrow(seasonMutDF2[seasonMutDF2$UMich=="UMich",]),
                           allMutatingSites,nchar(seasonMutDF2$HAH3muts[1]),
                           nrow(haploDF_wMich), nrow(haploDF_noMich),
                           paste(seasonMutDF2$highfreqMut[1], collapse= ","),
                           UMichEntropy, nonMichEntropy, allEntropy)
      }
    }
  }
}
resultsDF <- resultsDF[!is.na(resultsDF$subtype),]
write.csv(resultsDF, paste0(rootDir, "haplo_stats3/bigResultsDF_", Sys.Date(), ".csv"), row.names = F)

