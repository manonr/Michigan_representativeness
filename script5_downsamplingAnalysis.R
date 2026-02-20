# ============================================================================
# script5_downsamplingAnalysis.R
# Downsampling sensitivity analysis. Repeatedly subsamples Michigan sequences
# at various proportions, re-computes haplotype frequencies and detection
# stats across replicates, then aggregates results.
#
# Inputs: *_ALLmuts.rds
# Outputs: downsampled_haploStats_wReps_*.csv
# Dependencies: functions_pick_strains_functions.R
#
# KNOWN ISSUES:
#   - Calls entropyFct() which is commented out in functions_pick_strains_functions.R
#     (the v5 pipeline uses shannon_entropy instead -- this script may need updating)
#   - plot_functions.R is sourced but no plot functions are called
# ============================================================================

# version of script 1,  with downsampling reps

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
library(MMWRweek)


rootDir <- paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/")
setwd(rootDir)
dir.create("reps/")

source("/Users/manonragonnet/OneDrive/Github/Michigan_rep/R_Michigan/functions_pick_strains_functions.R")
source("/Users/manonragonnet/OneDrive/Github/Michigan_rep/R_Michigan/plot_functions.R")



subtypes <- c("H1", "H3")
thresholds <- c( 0.1)#0.05, 0.1, 0.2)
mySeasons <- c(2021, 2022, 2023, 2024, paste(2021:2024, collapse="-"))


UMich_ds <- c(0.05, 0.1, 0.25, 0.5, 0.75) # proportion of sequences you want to keep
reps <- 10

resultsDF_reps <- data.frame(subtype = "H1", threshold = 0.05, mySeasons = 2021,ds = 0.5, rep=1,
                             nSeq = 0, nMichSeq= 0,nVarSites= 0, nhighFreqSites= 0, 
                             nHaplo_wMich = 0, nHaplo_noMich = 0, highfreqMutSites="none",
                             UMichEntropy=0, nonMichEntropy = 0, allEntropy=0)
i <- 0
for (subtype in subtypes){
  #read once
  bigMutDF <- readRDS(paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/", subtype,"_ALLmuts.rds"))
  cat("[script5] ", subtype, " | loaded RDS | total:", nrow(bigMutDF),
      "| UMich:", sum(bigMutDF$UMich == "UMich"), "| notUOM:", sum(bigMutDF$UMich == "notUOM"), "\n")

  # do the ds here, across seasons, so that we can also check across years when haplotypes are found
  for (ds in UMich_ds){
    
    writeDir <- paste0(rootDir, "/reps/michDownsample", ds, "/")
    if(!dir.exists(writeDir)){
      dir.create(writeDir, recursive = T)
    }
    
    setwd(writeDir)
    
    
    for (rep in 1:reps){
      print(rep)
      
      bigMutDF_ds <- drop_michProp(bigMutDF, ds)
      cat("[script5] ", subtype, " | ds", ds, "rep", rep,
          "| after downsample | total:", nrow(bigMutDF_ds),
          "| UMich:", sum(bigMutDF_ds$UMich == "UMich"),
          "| notUOM:", sum(bigMutDF_ds$UMich == "notUOM"), "\n")

      write.csv(bigMutDF_ds, paste0(subtype, "_ds_", ds, "_rep_", rep, ".csv"), row.names = F) # for using later for dectecting haplotypes in later seasons
      
      for (threshold in thresholds){
        for (mySeason in mySeasons){
          if(length(strsplit(as.character(mySeason), "-")[[1]])>1){
            mySeason <- base::strsplit(as.character(mySeason), "-")[[1]]
          }
          
          runName <- paste0(subtype,"_ds", ds, "_rep", rep, "_t", threshold, "_", paste0(mySeason, collapse="."))
          print(runName)
          
          seasonMutDF <- bigMutDF_ds %>% filter(season %in% mySeason)
          allMutatingSites <- length((getMutList(seasonMutDF, "aapos", myCutoff=1/nrow(seasonMutDF))))
          seasonMutDF2 <- create_mut_df(seasonMutDF,threshold = threshold, mySeason = mySeason)
          
          if(!is.null(seasonMutDF2)){
            cat("[script5] ", runName,
                "| total:", nrow(seasonMutDF2),
                "| UMich:", sum(seasonMutDF2$UMich == "UMich"),
                "| notUOM:", sum(seasonMutDF2$UMich == "notUOM"), "\n")
            i<- i+1
            haploDF_noMich <- haploOrder(seasonMutDF2, inclMich = F)
            haploDF_wMich <- haploOrder(seasonMutDF2, inclMich = T)
            
            write.csv(haploDF_noMich, paste0(writeDir, "haploDF_noMich_", runName,".csv"), row.names = F)
            write.csv(haploDF_wMich, paste0(writeDir, "haploDF_wMich_", runName,".csv"), row.names = F)
            
            #entropy
            if(length(seasonMutDF2$HAH3muts[seasonMutDF2$UMich=="UMich"])>0){
              UMichHaplos <- as.matrix(AAStringSet(seasonMutDF2$HAH3muts[seasonMutDF2$UMich=="UMich"]))
              UMichEntropy <- round(mean(apply(UMichHaplos, 2, entropyFct)),3)
            } else {
              UMichHaplos <- 0
              UMichEntropy <- 0
            }
            
            
            otherHaplos <- as.matrix(AAStringSet(seasonMutDF2$HAH3muts[seasonMutDF2$UMich=="notUOM"]))
            allHaplos <- as.matrix(AAStringSet(seasonMutDF2$HAH3muts))
            
            nonMichEntropy <- round(mean(apply(otherHaplos, 2, entropyFct)),3)
            allEntropy <- round(mean(apply(allHaplos, 2, entropyFct)),3)
            
            # stats and results table
            resultsDF_reps[i,] <- c(subtype, threshold, paste(mySeason, collapse = "."),
                                    ds , rep,
                                    nrow(seasonMutDF2), nrow(seasonMutDF2[seasonMutDF2$UMich=="UMich",]),
                                    allMutatingSites,nchar(seasonMutDF2$HAH3muts[1]),
                                    nrow(haploDF_wMich), nrow(haploDF_noMich),
                                    paste(seasonMutDF2$highfreqMut[1], collapse= ","),
                                    UMichEntropy, nonMichEntropy, allEntropy)
          }
        }
      }
    }
  }
}
write.csv(resultsDF_reps, paste0(rootDir, "/reps/bigResultsDF_wReps_", Sys.Date(), ".csv"), row.names = F)





###############
# 2. detect haplos later -  from script 2
##############


rootDir <- paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/reps")
setwd(rootDir)


resultDFs <- dir()[grep("bigResultsDF_wReps", dir())]
bigResults <- read.csv(paste0(resultDFs[which.max(as.Date(sub("bigResultsDF_wReps_(\\d{4}-\\d{2}-\\d{2})\\.csv", "\\1", basename(resultDFs))))]))
bigResults <- bigResults[!is.na(bigResults$subtype),]

subtypes <- unique(bigResults$subtype)
thresholds <- unique(bigResults$threshold)
UMich_ds <- unique(bigResults$ds)
reps <- max(bigResults$rep)
mySeasons <- 2021:2024


for (subtype in subtypes){
  #bigMutDF <- readRDS(paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/", subtype,"_ALLmuts.rds"))
  for (ds in UMich_ds){
    setwd(paste0(rootDir, "/michDownsample",ds))
    dir.create("haplos_perSeason/")
    for (rep in 1:reps){
      #grab the correct down samples data
      
      ds_mutDF <- read.csv(paste0(subtype, "_ds_", ds, "_rep_", rep, ".csv"))
      
      
      for (threshold in thresholds){
        for (mySeason in mySeasons){
          runName <- paste0(subtype,"_ds", ds, "_rep", rep, "_t", threshold, "_", paste0(mySeason, collapse="."))
          print(runName)
          
          highFreqSites <- bigResults$highfreqMutSites[bigResults$subtype==subtype & 
                                                         bigResults$threshold==threshold & bigResults$mySeasons==mySeason
                                                       & bigResults$ds==ds & bigResults$rep==rep]
          if(length(highFreqSites)>0){
            highFreqSites <- as.numeric(strsplit(highFreqSites, " ")[[1]])
            mutDF <- create_mut_df_fromList(ds_mutDF, listMut = highFreqSites) # subset every sequence to the set of hgih freq sites
            haplos <- mutDF[,c("label","date", "season", "UMich", "HAH3muts", "location", "highfreqMut" )]
            write.csv(haplos, paste0("haplos_perSeason/",runName,"_haplos.csv"), row.names = F)
          }
        }
      }
    }
  }
}



###############
# 3. Haplotype stats - from script 3
##############



rootDir <- paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/reps/")
setwd(rootDir)



mySeasons <- unique(bigResults$mySeasons)


haploStatsDF <- data.frame(subtype =NA, threshold=0, season = NA, downsampleProp = 1, rep=1,
                           p90noMich=0,p50noMich=0,topXinMich_noMich = 0,ofTop50_nMissingInMich_noMich = 0,ofTop90_nMissingInMich_noMich = 0,
                           sumPropStrains_inMich_noMich = 0,sumPropStrains_missingMich_noMich = 0,
                           p90Mich=0, p50Mich=0, topXinMich_wMich = 0, ofTop50_nMissingInMich_wMich = 0, 
                           ofTop90_nMissingInMich_wMich = 0, sumPropStrains_inMich_wMich = 0, sumPropStrains_missingMich_wMich = 0)
i <- 0
for (subtype in subtypes){
  for (ds in UMich_ds){
    setwd(paste0(rootDir, "michDownsample", ds))
    
    for (rep in 1:reps){
      for (threshold in thresholds){
        for (mySeason in mySeasons){
          
          i <- i+1
          runName <- paste0(subtype,"_ds", ds, "_rep", rep, "_t", threshold, "_", paste0(mySeason, collapse="."))
          print(runName)
          
          haploDF_noMich <- read.csv(paste0("haploDF_noMich_", runName,".csv"))
          haploDF_wMich <- read.csv(paste0("haploDF_wMich_", runName,".csv"))
          
          #update with haplos that were detected at a later date
          # identify the correct file
          if(mySeason != "2021.2022.2023.2024"){
            
            fullHaploinfo_fileName <- paste0("haplos_perSeason/",runName , "_haplos.csv")
            fullHaploinfo <- read.csv(fullHaploinfo_fileName)
            
            haploDF_wMich2 <- fill_in_earlierHaplos(haploDF_wMich, fullHaploinfo, mySeason = mySeason)
            haploDF_noMich2 <- fill_in_earlierHaplos(haploDF_noMich, fullHaploinfo, mySeason = mySeason)
            write.csv(haploDF_noMich2, paste0(rootDir, "michDownsample", ds, "/haploDF_noMich_", runName,"_v2.csv"), row.names = F)
            write.csv(haploDF_wMich2, paste0(rootDir, "michDownsample", ds, "/haploDF_wMich_", runName,"_v2.csv"), row.names = F)
          } else {
            haploDF_wMich2 <- haploDF_wMich
            haploDF_noMich2 <- haploDF_noMich
          }
          
          # what would be more informative, below, than the number of syrains not in michigan, is the proportion that they represent
          # X% of strains have a Michigan haplo match, X% f srains do not
          
          haploStatsDF[i,] <-  c(subtype, threshold, paste(mySeason, collapse = "."),ds, rep,
                                 unlist(haplotype_stats(haploDF_noMich2)), unlist(haplotype_stats(haploDF_wMich2)))
        }
      }
    }
  }
}
write.csv(haploStatsDF, paste0(rootDir, "bigHaploStatsDF_wReps_", Sys.Date(), ".csv"))



############
# change in delays and prop detected
###########

#want to show X axis : proportion sampled
#y axis : number of NAs, delay in detection



rootDir <- "/Users/manonragonnet/Library/CloudStorage/OneDrive-Personal/Projects/LauringLab/Representativeness/haplo_stats3/reps"
setwd(rootDir)

myFiles <- list.files(pattern = "haploDF_noMich.*_v2", recursive = TRUE)
all_data <- map_df(myFiles, ~ read_csv(.x) %>% mutate(source_file = basename(.x)))




all_data$season <- all_data$USseason_dominant
all_data$subtype <- unlist(lapply(all_data$source_file, function(x) {paste0("H", strsplit(strsplit(x, "_H")[[1]][2], "_")[[1]][1])}))
all_data$downsample <- sub(".*ds([0-9.]+).*", "\\1", all_data$source_file)
all_data$threshold <- sub(".*t([0-9.]+).*", "\\1", all_data$source_file)
all_data$rep<- sub(".*rep([0-9.]+).*", "\\1", all_data$source_file)


all_data <- all_data %>% mutate(detected  = ifelse(!is.na(delay),1,0))

# drop h1 2021
dropLines1 <- which(all_data$subtype=="H1" & all_data$season=="2021")
dropLines2 <- which(is.na(all_data$firstUSDate))  # weird result here with H1 2021 0.005 to check # fixed, should be empty

all_data <- all_data[- c(dropLines1, dropLines2),]
all_data$WEEK <- MMWRweek(all_data$firstUSDate)$MMWRweek
table(all_data$firstUS_season, all_data$subtype, all_data$threshold, all_data$downsample)

write.csv(all_data, paste0(rootDir, "/downsampled_haploStats_wReps_", Sys.Date(), ".csv"), row.names = F)


#15 min so far



