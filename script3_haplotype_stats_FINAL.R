# ============================================================================
# script3_haplotype_stats_v2.R
# Step 3: Compute haplotype detection statistics. Updates first-detection
# dates across all seasons and calculates summary stats (p50, p90 coverage,
# missing haplotypes in Michigan).
#
# Inputs: haploDF_*.csv, haplos_perSeason/*.csv
# Outputs: haploDF_*_v2.csv, bigHaploStatsDF_*.csv
# Dependencies: functions_pick_strains_functions.R
#
# KNOWN ISSUES:
#   - source() path is R_Michigan/functions_pick_strains_functions.R (missing
#     Representativeness/ subdirectory) -- may point to a different copy
# ============================================================================

# now with version for reps/ without reps

source("/Users/manonragonnet/OneDrive/Github/Michigan_rep/R_Michigan/functions_pick_strains_functions.R")


# haplotype stats
rootDir <- paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/")
setwd(rootDir)


resultDFs <- dir()[grep("bigResultsDF", dir())]
bigResults <- read.csv(paste0(resultDFs[which.max(as.Date(sub("bigResultsDF_(\\d{4}-\\d{2}-\\d{2})\\.csv", "\\1", basename(resultDFs))))]))
bigResults <- bigResults[!is.na(bigResults$subtype),]

subtypes <- unique(bigResults$subtype)
thresholds <- unique(bigResults$threshold)
mySeasons <- unique(bigResults$mySeasons)


haploStatsDF <- data.frame(subtype ="H1", threshold=0, season = 2022, 
                           p90noMich=0,p50noMich=0,topXinMich_noMich = 0,ofTop50_nMissingInMich_noMich = 0,ofTop90_nMissingInMich_noMich = 0,
                           sumPropStrains_inMich_noMich = 0,sumPropStrains_missingMich_noMich = 0,
                           p90Mich=0, p50Mich=0, topXinMich_wMich = 0, ofTop50_nMissingInMich_wMich = 0, 
                           ofTop90_nMissingInMich_wMich = 0, sumPropStrains_inMich_wMich = 0, sumPropStrains_missingMich_wMich = 0)
i <- 0
for (subtype in subtypes){
  for (threshold in thresholds){
    for (mySeason in mySeasons){
      
      runName <- paste(subtype, threshold, paste(mySeason, collapse = "."), sep="_")
      print(runName)
      
      if(file.exists(paste0(rootDir, "haploDF_wMich_", runName,".csv"))){
        i <- i+1
        
        haploDF_wMich <- read.csv(paste0(rootDir, "haploDF_wMich_", runName,".csv"))
        haploDF_noMich <- read.csv(paste0(rootDir, "haploDF_noMich_", runName,".csv"))
        cat("[script3] ", runName,
            "| wMich haplos:", nrow(haploDF_wMich), "total seqs:", sum(haploDF_wMich$Count),
            "| noMich haplos:", nrow(haploDF_noMich), "total seqs:", sum(haploDF_noMich$Count), "\n")

        #update with haplos that were detected EARLIER OR at a later date
        # identify the correct file
        if(mySeason != "2021.2022.2023.2024"){
          
          fullHaploinfo_fileName <- paste0("haplos_perSeason/", subtype, "_", threshold, "_", mySeason, "_haplos.csv")
          fullHaploinfo <- read.csv(fullHaploinfo_fileName)
          
          
          haploDF_wMich2 <- fill_in_earlierHaplos(haploDF_wMich, fullHaploinfo, mySeason = mySeason)
          haploDF_noMich2 <- fill_in_earlierHaplos(haploDF_noMich, fullHaploinfo, mySeason = mySeason)
          write.csv(haploDF_wMich2, paste0(rootDir, "haploDF_wMich_", runName,"_v2.csv"), row.names = F)
          write.csv(haploDF_noMich2, paste0(rootDir, "haploDF_noMich_", runName,"_v2.csv"), row.names = F)
        } else {
          haploDF_wMich2 <- haploDF_wMich
          haploDF_noMich2 <- haploDF_noMich
        }
        
        # what would be more informative, below, than the number of syrains not in michigan, is the proportion that they represent
        # X% of strains have a Michigan haplo match, X% f srains do not
        
        haploStatsDF[i,] <-  c(subtype, threshold, paste(mySeason, collapse = "."),
                               unlist(haplotype_stats(haploDF_noMich2)), unlist(haplotype_stats(haploDF_wMich2)))
        
      }
    }
  }
}
write.csv(haploStatsDF, paste0(rootDir, "bigHaploStatsDF_", Sys.Date(), ".csv"))

