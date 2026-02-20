# ============================================================================
# script7_otherStateProxies.R
# State-level proxy analysis. Extends the haplotype representativeness
# analysis to all US states and compares each state's detection performance
# against Michigan.
#
# Inputs: *_ALLmuts.rds
# Outputs: bigResultsDF_stateProxyAnalysis_*.csv,
#          haploStats_stateProxyAnalysis_*.csv
# Dependencies: functions_pick_strains_functions.R, state_analysis_functions.R
#
# KNOWN ISSUES:
#   - read_one() function is defined twice (identical copies)
#   - plot_functions.R is sourced but no plot functions are called
#   - Commented-out entropy code blocks (lines ~82-96)
# ============================================================================

# version of script 1,  with state analysis

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
library(purrr)
library(pracma) 


rootDir <- paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/")
setwd(rootDir)

source("/Users/manonragonnet/OneDrive/Github/Michigan_rep/R_Michigan/functions_pick_strains_functions.R")
source("/Users/manonragonnet/OneDrive/Github/Michigan_rep/R_Michigan/plot_functions.R")
source("/Users/manonragonnet/OneDrive/Github/Michigan_rep/R_Michigan/state_analysis_functions.R")



subtypes <- c("H1", "H3")
thresholds <- c( 0.1)#, 0.05, 0.1, 0.2)
mySeasons <- c(2021, 2022, 2023, 2024, paste(2021:2024, collapse="-"))

writeDir <- paste0(rootDir, "stateProxyAnalysis/")
if(!dir.exists(writeDir)){
  dir.create(writeDir)
}




resultsDF_epi <- data.frame(subtype = "H1", threshold = 0.05, mySeasons = 2021,state = "US",
                            nSeq = 0, nStateSeq= 0,nVarSites= 0, nhighFreqSites= 0, 
                            nHaplo_noState = 0, highfreqMutSites="none")
# UMichEntropy=0, nonMichEntropy = 0, allEntropy=0)
i <- 0


for (subtype in subtypes){
  #read once
  bigMutDF <- readRDS(paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/", subtype,"_ALLmuts.rds"))
  cat("[script7] ", subtype, " | loaded RDS | total:", nrow(bigMutDF),
      "| UMich:", sum(bigMutDF$UMich == "UMich"), "| notUOM:", sum(bigMutDF$UMich == "notUOM"), "\n")

  setwd(writeDir)
  
  for (threshold in thresholds){
    for (mySeason in mySeasons){
      if(length(strsplit(as.character(mySeason), "-")[[1]])>1){
        mySeason <- base::strsplit(as.character(mySeason), "-")[[1]]
      }
      
      
      seasonMutDF <- bigMutDF %>% filter(season %in% mySeason) %>% dplyr::rename(state=location)
      allMutatingSites <- length((getMutList(seasonMutDF, "aapos", myCutoff=1/nrow(seasonMutDF))))
      seasonMutDF2 <- create_mut_df(seasonMutDF,threshold=threshold, mySeason=mySeason) # only sites with some variatiojn are in seaosnMutDF, so need to subset to waht is available there
      
      if(!is.null(seasonMutDF2)){
        cat("[script7] ", subtype, " | threshold", threshold, "| season", paste(mySeason, collapse="."),
            "| total:", nrow(seasonMutDF2),
            "| UMich:", sum(seasonMutDF2$UMich == "UMich"),
            "| notUOM:", sum(seasonMutDF2$UMich == "notUOM"),
            "| nStates:", length(unique(seasonMutDF2$state)), "\n")

        bigStates <- unique(seasonMutDF2$state) #names(table(seasonMutDF2$state)[which(table(seasonMutDF2$state)>=20)])
        
        for (state in bigStates){
          
          runName <- paste0(subtype,"_", state, "_stateProxyAnalysis_t", threshold, "_", paste0(mySeason, collapse="."))
          print(runName)
          
          i<- i+1
          haploDF_state <- haploOrder_pickState(seasonMutDF2, state = state)
          write.csv(haploDF_state, paste0(writeDir, "haploDF_",  runName,".csv"), row.names = F)
          
          ##entropy
          #if(length(seasonMutDF2$HAH3muts[seasonMutDF2$UMich=="UMich"])>0){
          # UMichHaplos <- as.matrix(AAStringSet(seasonMutDF2$HAH3muts[seasonMutDF2$UMich=="UMich"]))
          #  UMichEntropy <- round(mean(apply(UMichHaplos, 2, entropyFct)),3)
          #} else {
          #  UMichHaplos <- 0
          #  UMichEntropy <- 0
          #}
          
          
          #otherHaplos <- as.matrix(AAStringSet(seasonMutDF2$HAH3muts[seasonMutDF2$UMich=="notUOM"]))
          allHaplos <- as.matrix(AAStringSet(seasonMutDF2$HAH3muts))
          
          #nonMichEntropy <- round(mean(apply(otherHaplos, 2, entropyFct)),3)
          #allEntropy <- round(mean(apply(allHaplos, 2, entropyFct)),3)
          
          # stats and results table
          resultsDF_epi[i,] <- c(subtype, threshold, paste(mySeason, collapse = "."),state,
                                 nrow(seasonMutDF2), nrow(seasonMutDF2[seasonMutDF2$state==state,]),
                                 allMutatingSites,nchar(seasonMutDF2$HAH3muts[1]),
                                 nrow(haploDF_state),
                                 paste(seasonMutDF2$highfreqMut[1], collapse= ","))
          #UMichEntropy, nonMichEntropy, allEntropy)
        }
      }
    }
  }
}

write.csv(resultsDF_epi, paste0(writeDir, "bigResultsDF_stateProxyAnalysis_", Sys.Date(), ".csv"), row.names = F)





###############
# 2. detect haplos later -  from script 2
##############


rootDir <- paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/stateProxyAnalysis/")
setwd(rootDir)


resultDFs <- dir()[grep("bigResultsDF_stateProxyAnalysis", dir())]
bigResults <- read.csv(paste0(resultDFs[which.max(as.Date(sub("bigResultsDF_stateProxyAnalysis_(\\d{4}-\\d{2}-\\d{2})\\.csv", "\\1", basename(resultDFs))))]))
bigResults <- bigResults[!is.na(bigResults$subtype),]

subtypes <- unique(bigResults$subtype)
thresholds <- unique(bigResults$threshold)
mySeasons <- unique(bigResults$mySeasons)#2021:2024


for (subtype in subtypes){
  bigMutDF <- readRDS(paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/", subtype,"_ALLmuts.rds"))
  setwd(paste0(rootDir))
  dir.create("haplos_perSeason/")
  #grab the correct down samples data
  
  for (threshold in thresholds){
    for (mySeason in mySeasons){
      runName_short <- paste0(subtype,"_stateProxyAnalysis_t", threshold, "_", paste0(mySeason, collapse="."))
      print(runName_short)
      
      highFreqSites <- bigResults$highfreqMutSites[bigResults$subtype==subtype & 
                                                     bigResults$threshold==threshold & bigResults$mySeasons==mySeason]
      if(length(highFreqSites)>0){
        highFreqSites <- as.numeric(strsplit(highFreqSites, " ")[[1]])
        mutDF <- create_mut_df_fromList(bigMutDF, listMut = highFreqSites) # subset every sequence to the set of hgih freq sites
        haplos <- mutDF[,c("label","date", "season", "UMich", "HAH3muts", "location", "highfreqMut" )]
        write.csv(haplos, paste0("haplos_perSeason/",runName_short,"_haplos.csv"), row.names = F)
      }
    }
  }
}




###############
# 3. Haplotype stats - from script 3
##############



rootDir <- paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/stateProxyAnalysis/")

setwd(rootDir)

mySeasons <- unique(bigResults$mySeasons)


haploStatsDF <- data.frame(subtype =NA, threshold=0, season = NA, state=NA,
                           p90=0,p50=0,topXinState = 0,ofTop50_nMissingInState = 0,ofTop90_nMissingInState = 0,
                           sumPropStrains_inState = 0,sumPropStrains_missingState = 0)
i <- 0
for (thisSubtype in subtypes){
  
  setwd(rootDir)
  
  for (thisThreshold in thresholds){
    for (thisSeason in mySeasons){
      
      runName_short <- paste0(thisSubtype,"_stateProxyAnalysis_t", thisThreshold, "_", paste0(thisSeason, collapse="."))
      print(runName_short)
      #subset df to see which states were processed
      states <- bigResults %>%  dplyr::filter(subtype == thisSubtype, threshold == thisThreshold, mySeasons == thisSeason) %>%pull(state)
      states <- unique(states)
      
      for (state in states){
        
        
        i <- i+1
        runName <- paste0(thisSubtype,"_", state, "_stateProxyAnalysis_t", thisThreshold, "_", paste0(thisSeason, collapse="."))
        print(runName)
        
        haploDF <- read.csv(paste0("haploDF_", runName,".csv"))
        
        
        #update with haplos that were detected at a later date
        # identify the correct file
        # if(mySeason != "2021.2022.2023.2024"){
        
        fullHaploinfo_fileName <- paste0("haplos_perSeason/",runName_short , "_haplos.csv")
        fullHaploinfo <- read.csv(fullHaploinfo_fileName)
        names(fullHaploinfo)[ names(fullHaploinfo)=="location"] <- "state"
        
        haploDF_v2 <- fill_in_earlierHaplos_allStates(haploDF, fullHaploinfo, mySeason=thisSeason, mystate = state)
        
        
        write.csv(haploDF_v2, paste0(rootDir,  "/haploDF_", runName,"_v2.csv"), row.names = F)
        
        # what would be more informative, below, than the number of syrains not in michigan, is the proportion that they represent
        # X% of strains have a Michigan haplo match, X% f srains do not
        
        haploStatsDF[i,] <-  c(thisSubtype, thisThreshold, paste(thisSeason, collapse = "."),state,
                               unlist(haplotype_stats_otherStates(haploDF_v2)))      }
    }
  }
}

write.csv(haploStatsDF, paste0(rootDir, "bigHaploStatsDF_stateProxyAnalysis_", Sys.Date(), ".csv"), row.names = F)



############
# change in delays and prop detected
###########

#want to show X axis : proportion sampled
#y axis : number of NAs, delay in detection



rootDir <- "/Users/manonragonnet/Library/CloudStorage/OneDrive-Personal/Projects/LauringLab/Representativeness/haplo_stats3/stateProxyAnalysis/"
setwd(rootDir)

myFiles <- list.files(pattern = "haploDF.*_v2", recursive = TRUE)
#myFiles <- myFiles[grep("Michigan", myFiles)]

read_one <- function(p) {
  read_csv(p, col_types = cols(USseason_dominant = col_character())) |>
    mutate(source_file = basename(p))
}

all_data <- map_dfr(myFiles, read_one)


all_data$season <- all_data$USseason_dominant
all_data$subtype <- unlist(lapply(all_data$source_file, function(x) {paste0("H", strsplit(strsplit(x, "_H")[[1]][2], "_")[[1]][1])}))
all_data$threshold <- sub(".*t([0-9.]+).*", "\\1", all_data$source_file)



all_data <- all_data %>% mutate(detected  = ifelse(!is.na(delay),1,0))

# drop h1 2021 and weird
table(all_data$USseason_dominant, all_data$subtype, all_data$threshold)
dropLines1 <- which(all_data$subtype=="H1" & all_data$USseason_dominant=="2021")
dropLines2 <- which(is.na(all_data$firstUSDate))  # weird result here with H1 2021 0.005 to check # fixed, should be empty

all_data <- all_data[- c(dropLines1, dropLines2),]
all_data$WEEK <- MMWRweek(all_data$firstUSDate)$MMWRweek
table(all_data$firstUS_season, all_data$subtype, all_data$threshold)

write.csv(all_data, paste0(rootDir, "/haploStats_stateProxyAnalysis_", Sys.Date(), ".csv"), row.names = F)



# stat comparison btw each state and michigan

reps <- 10

rootDir <- paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/stateProxyAnalysis/")
setwd(rootDir)

writeDir <- paste0(rootDir, "reps")
if(!dir.exists(writeDir)){
  dir.create(writeDir)
}


resultDFs <- dir()[grep("bigResultsDF_stateProxyAnalysis", dir())]
bigResults <- read.csv(paste0(resultDFs[which.max(as.Date(sub("bigResultsDF_stateProxyAnalysis_(\\d{4}-\\d{2}-\\d{2})\\.csv", "\\1", basename(resultDFs))))]))
bigResults
bigResults <- bigResults[complete.cases(bigResults),]#bigResults[!is.na(bigResults$subtype),]


resultsDF_reps <- data.frame(subtype = "H1", state = "Alabama", threshold = 0.05, mySeasons = 2021, rep=1, nStateSeq=0,
                             nSeq = 0, nMichSeq= 0,nVarSites= 0, nhighFreqSites= 0, 
                             nHaplo_noMich = 0, highfreqMutSites="none")

i <- 0
for (subtype in unique(bigResults$subtype)){
  
  bigMutDF <- readRDS(paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/",  subtype,"_ALLmuts.rds"))
  
  for (threshold in unique(bigResults$threshold)){
    for (state in unique(bigResults$state)){
      
      for (rep in 1:reps){
        
        # have to sample across all seasons, with numbers that match each season independently
        samplingFreq <- bigResults[bigResults$subtype==subtype & bigResults$threshold==threshold & bigResults$state==state,c("mySeasons", "nStateSeq")]
        print(samplingFreq)
        samplingFreq <- samplingFreq[samplingFreq$mySeasons %in% 2021:2024,]
        bigMutDF_ds <- keepMich_n(mutDF=bigMutDF, mySamplingFreq=samplingFreq)
        cat("[script7-reps] ", subtype, " | state", state, "| rep", rep,
            "| after downsample | total:", nrow(bigMutDF_ds),
            "| UMich:", sum(bigMutDF_ds$UMich == "UMich"),
            "| notUOM:", sum(bigMutDF_ds$UMich == "notUOM"), "\n")
        write.csv(bigMutDF_ds, paste0("reps/",subtype, "_Mich_ds_", state, "_rep", rep, "_t", threshold,".csv"), row.names = F)
        
        for (mySeason in unique(bigResults$mySeasons)){
          
          if(length(strsplit(as.character(mySeason), "\\.")[[1]])>1){
            mySeason <- base::strsplit(as.character(mySeason), "\\.")[[1]]
          }
          
          runName <- paste0(subtype,"_", state,
                            "_rep", rep, "_t", threshold ,"_", paste0(mySeason, collapse="."))
          print(runName)
          
          stateSeq <- samplingFreq$nStateSeq[samplingFreq$mySeasons==paste0(mySeason, collapse=".")]
          print("Note in current form, I am not running stats on multi season data, only season by season")
          
          if(length(stateSeq>0)){
            
            
            seasonMutDF <- bigMutDF_ds %>% filter(season %in% mySeason)
            allMutatingSites <- length((getMutList(seasonMutDF, "aapos", myCutoff=1/nrow(seasonMutDF))))
            seasonMutDF2 <- create_mut_df(seasonMutDF,threshold = threshold, mySeason = mySeason)
            
            if(!is.null(seasonMutDF2)){
              print("c")
              i<- i+1
              haploDF_noMich <- haploOrder(seasonMutDF2, inclMich = F)
              write.csv(haploDF_noMich, paste0(writeDir, "/haploDF_noMich_", runName,".csv"), row.names = F)
              resultsDF_reps[i,] <- c(subtype, state, threshold, paste(mySeason, collapse = "."),
                                      rep,
                                      stateSeq , nrow(seasonMutDF2), nrow(seasonMutDF2[seasonMutDF2$UMich=="UMich",]),
                                      allMutatingSites,nchar(seasonMutDF2$HAH3muts[1]),
                                      nrow(haploDF_noMich),
                                      paste(seasonMutDF2$highfreqMut[1], collapse= ","))
            }
            
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

setwd(paste0(rootDir, "/reps"))


resultDFs <- dir()[grep("bigResultsDF_wReps", dir())]
bigResults <- read.csv(paste0(resultDFs[which.max(as.Date(sub("bigResultsDF_wReps_(\\d{4}-\\d{2}-\\d{2})\\.csv", "\\1", basename(resultDFs))))]))
bigResults <- bigResults[!is.na(bigResults$subtype),]

subtypes <- unique(bigResults$subtype)
thresholds <- unique(bigResults$threshold)
states <- "Michigan"#unique(bigResults$state)
reps <- max(bigResults$rep)
mySeasons <- 2021:2024


for (subtype in subtypes){
  #bigMutDF <- readRDS(paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/", subtype,"_ALLmuts.rds"))
  dir.create("haplos_perSeason/")
  for (state in states){
    for (threshold in thresholds){
      for (rep in 1:reps){
        #grab the correct down samples data
        
        ds_mutDF <- read.csv(paste0(subtype, "_Mich_ds_", state, "_rep", rep, "_t", threshold,".csv"))
        
        
        
        for (mySeason in mySeasons){
          runName <- paste0(subtype,"_", state, "_rep", rep, "_t", threshold, "_", paste0(mySeason, collapse="."))
          print(runName)
          
          highFreqSites <- bigResults$highfreqMutSites[bigResults$subtype==subtype & 
                                                         bigResults$threshold==threshold & bigResults$mySeasons==mySeason
                                                       & bigResults$state==state & bigResults$rep==rep]
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



setwd(paste0(rootDir, "/reps/"))

mySeasons <- unique(bigResults$mySeasons)


haploStatsDF <- data.frame(subtype =NA, threshold=0, season = NA, state = NA, rep=1,
                           p90noMich=0,p50noMich=0,topXinMich_noMich = 0,ofTop50_nMissingInMich_noMich = 0,ofTop90_nMissingInMich_noMich = 0,
                           sumPropStrains_inMich_noMich = 0,sumPropStrains_missingMich_noMich = 0)
i <- 0
for (subtype in subtypes){
  for (state in states){
    for (threshold in thresholds){
      for (rep in 1:reps){
        
        for (mySeason in mySeasons){
          
          i <- i+1
          runName <- paste0(subtype,"_", state, "_rep", rep, "_t", threshold, "_", paste0(mySeason, collapse="."))
          print(runName)
          
          if(file.exists(paste0("haploDF_noMich_", runName,".csv"))){
            
            
            haploDF_noMich <- read.csv(paste0("haploDF_noMich_", runName,".csv"))
            
            
            #update with haplos that were detected at a later date
            # identify the correct file
            if(mySeason != "2021.2022.2023.2024"){
              
              fullHaploinfo_fileName <- paste0("haplos_perSeason/",runName , "_haplos.csv")
              fullHaploinfo <- read.csv(fullHaploinfo_fileName)
              
              haploDF_noMich2 <- fill_in_earlierHaplos(haploDF_noMich, fullHaploinfo, mySeason = mySeason)
              
            } else {
              haploDF_noMich2 <- haploDF_noMich
            }
            
            write.csv(haploDF_noMich2, paste0(rootDir, "/reps/haploDF_noMich_", runName,"_v2.csv"), row.names = F)
            
            
            # what would be more informative, below, than the number of syrains not in michigan, is the proportion that they represent
            # X% of strains have a Michigan haplo match, X% f srains do not
            
            haploStatsDF[i,] <-  c(subtype, threshold, paste(mySeason, collapse = "."),state, rep,
                                   unlist(haplotype_stats(haploDF_noMich2)))
          }
        }
      }
    }
  }
}
write.csv(haploStatsDF, paste0(rootDir, "/reps/bigHaploStatsDF_wReps_", Sys.Date(), ".csv"))

# combine all data for the reps analysis

setwd(paste0(rootDir, "reps/"))

myFiles <- list.files(pattern = "haploDF.*_v2", recursive = TRUE)
#myFiles <-myFiles[grep("Michigan", myFiles)]

read_one <- function(p) {
  read_csv(p, col_types = cols(USseason_dominant = col_character())) |>
    mutate(source_file = basename(p))
}

all_data2 <- map_dfr(myFiles, read_one)


all_data2$season <- all_data2$USseason_dominant
all_data2$subtype <- unlist(lapply(all_data2$source_file, function(x) {paste0("H", strsplit(strsplit(x, "_H")[[1]][2], "_")[[1]][1])}))
all_data2$threshold <- sub(".*t([0-9.]+).*", "\\1", all_data2$source_file)
all_data2$rep <- sub(".*rep([0-9.]+).*", "\\1", all_data2$source_file)
all_data2$thisState <- unlist(lapply(all_data2$source_file, function(x) {strsplit(x, "_")[[1]][4]}))
all_data2 <- all_data2 %>% mutate(detected  = ifelse(!is.na(delay),1,0))

write.csv(all_data2, paste0(rootDir, "/reps/haploStats_stateProxyAnalysis_MichReps_", Sys.Date(), ".csv"), row.names = F)


