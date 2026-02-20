# ============================================================================
# script2_detecting_haplos_later_v2.R
# Step 2: Track haplotypes across seasons. Applies mutation site definitions
# from Step 1 to ALL seasons, enabling detection of when each haplotype
# first appeared.
#
# Inputs: bigResultsDF_*.csv, *_ALLmuts.rds
# Outputs: haplos_perSeason/*.csv
# Dependencies: functions_pick_strains_functions.R
# ============================================================================


# detecting delays requires running all sequences with all haplotype definitions, then updating the delay table

source("/Users/manonragonnet/OneDrive/Github/Michigan_rep/R_Michigan/Representativeness/functions_pick_strains_functions.R")


rootDir <- paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/")
setwd(rootDir)
dir.create("haplos_perSeason/")

resultDFs <- dir()[grep("bigResultsDF", dir())]
bigResults <- read.csv(paste0(resultDFs[which.max(as.Date(sub("bigResultsDF_(\\d{4}-\\d{2}-\\d{2})\\.csv", "\\1", basename(resultDFs))))]))
bigResults <- bigResults[!is.na(bigResults$subtype),]

subtypes <- unique(bigResults$subtype)
thresholds <- unique(bigResults$threshold)
mySeasons <- 2021:2024


for (subtype in subtypes){
  bigMutDF <- readRDS(paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/", subtype,"_ALLmuts.rds"))
  cat("[script2] ", subtype, " | loaded RDS | total:", nrow(bigMutDF),
      "| UMich:", sum(bigMutDF$UMich == "UMich"), "| notUOM:", sum(bigMutDF$UMich == "notUOM"), "\n")
  for (threshold in thresholds){
    for (mySeason in mySeasons){
      print(paste(subtype, threshold, mySeason))
      highFreqSites <- bigResults$highfreqMutSites[bigResults$subtype==subtype & bigResults$threshold==threshold & bigResults$mySeasons==mySeason]
      if(length(highFreqSites)>0){
        highFreqSites <- as.numeric(strsplit(highFreqSites, " ")[[1]])
        mutDF <- create_mut_df_fromList(bigMutDF, listMut = highFreqSites)%>% dplyr::rename(state=location) # subset every sequence to the set of hgih freq sites
        cat("[script2] ", subtype, " | threshold", threshold, "| season", mySeason,
            "| mutDF total:", nrow(mutDF),
            "| UMich:", sum(mutDF$UMich == "UMich"),
            "| notUOM:", sum(mutDF$UMich == "notUOM"), "\n")
        runName <- paste(subtype, threshold, paste(mySeason, collapse = "."), sep="_")
        haplos <- mutDF[,c("label","date", "season", "UMich", "HAH3muts", "state", "highfreqMut" )]
        write.csv(haplos, paste("haplos_perSeason/",subtype, "_",threshold, "_",paste(mySeason, collapse = "."),"_haplos.csv", sep=""), row.names = F)
      }
    }
  }
}
  
  