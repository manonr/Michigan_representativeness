# ============================================================================
# script6_epitopeAnalysis.R
# Epitope-restricted haplotype analysis. Same pipeline as the main analysis
# but restricted to structural epitope sites on HA.
#
# Inputs: *_AA_structuralEpitopes.fas, *_ALLmuts.rds
# Outputs: bigResultsDF_wEpitopes_*.csv, haploStats_wEpitopes_*.csv
# Dependencies: functions_pick_strains_functions.R
#
# KNOWN ISSUES:
#   - [FIXED] !exists(writeDir) → !dir.exists(writeDir)
#   - [FIXED] Duplicate threshold 0.1 removed from vector
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


rootDir <- paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/")
setwd(rootDir)

source("/Users/manonragonnet/OneDrive/Github/Michigan_rep/R_Michigan/functions_pick_strains_functions.R")
source("/Users/manonragonnet/OneDrive/Github/Michigan_rep/R_Michigan/plot_functions.R")



subtypes <- c("H1", "H3")
thresholds <- c(0.05, 0.1, 0.2) #claude-fixed-error: was c(0.1, 0.05, 0.1, 0.2) -- 0.1 was duplicated
mySeasons <- c(2021, 2022, 2023, 2024, paste(2021:2024, collapse="-"))


epitopeSiteList <- function(subtype){
  setwd("/Users/manonragonnet/Library/CloudStorage/OneDrive-Personal/Projects/CobeyLab/HIVE_Michigan/library")
  #read in fasta
  seqFileName <- paste0(subtype, "_AA_structuralEpitopes.fas")
  aa_seq <- read.FASTA(seqFileName, type="AA")
  chars <- sapply(aa_seq$Epitopes, function(x) rawToChar(as.raw(strtoi(x, 16))))
  positions <- which(chars != "-")
}



# no difference in the list between years so maybe analysis does not need ot be done seaosn by season


resultsDF_epi <- data.frame(subtype = "H1", threshold = 0.05, mySeasons = 2021,
                            nSeq = 0, nMichSeq= 0,nVarSites= 0, nhighFreqSites= 0, 
                            nHaplo_wMich = 0, nHaplo_noMich = 0, highfreqMutSites="none",
                            UMichEntropy=0, nonMichEntropy = 0, allEntropy=0)
i <- 0
for (subtype in subtypes){
  #read once
  bigMutDF <- readRDS(paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/", subtype,"_ALLmuts.rds"))
  cat("[script6] ", subtype, " | loaded RDS | total:", nrow(bigMutDF),
      "| UMich:", sum(bigMutDF$UMich == "UMich"), "| notUOM:", sum(bigMutDF$UMich == "notUOM"), "\n")

  writeDir <- paste0(rootDir, "haplo_stats3/epitopeAnalysis/")
  if(!dir.exists(writeDir)){ #claude-fixed-error: was !exists() which checks R variables, not filesystem
    dir.create(writeDir)
  }
  
  setwd(writeDir)
  
  for (threshold in thresholds){
    for (mySeason in mySeasons){
      if(length(strsplit(as.character(mySeason), "-")[[1]])>1){
        mySeason <- base::strsplit(as.character(mySeason), "-")[[1]]
      }
      
      runName <- paste0(subtype,"_epitopeAnalysis_t", threshold, "_", paste0(mySeason, collapse="."))
      print(runName)
      
      seasonMutDF <- bigMutDF %>% filter(season %in% mySeason)
      epiSiteList<- epitopeSiteList(subtype)
      
      seasonMutDF2 <- create_mut_df_fromList(seasonMutDF,epiSiteList[paste0("aapos", epiSiteList) %in% names(seasonMutDF)]) # only sites with some variatiojn are in seaosnMutDF, so need to subset to waht is available there
      
      if(!is.null(seasonMutDF2)){
        cat("[script6] ", runName,
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
        resultsDF_epi[i,] <- c(subtype, threshold, paste(mySeason, collapse = "."),
                               nrow(seasonMutDF2), nrow(seasonMutDF2[seasonMutDF2$UMich=="UMich",]),
                               length(epiSiteList),nchar(seasonMutDF2$HAH3muts[1]),
                               nrow(haploDF_wMich), nrow(haploDF_noMich),
                               paste(seasonMutDF2$highfreqMut[1], collapse= ","),
                               UMichEntropy, nonMichEntropy, allEntropy)
      }
    }
  }
}

write.csv(resultsDF_epi, paste0(rootDir, "haplo_stats3/epitopeAnalysis/bigResultsDF_wEpitopes_", Sys.Date(), ".csv"), row.names = F)





###############
# 2. detect haplos later -  from script 2
##############


rootDir <- paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/epitopeAnalysis/")
setwd(rootDir)


resultDFs <- dir()[grep("bigResultsDF_wEpitopes", dir())]
bigResults <- read.csv(paste0(resultDFs[which.max(as.Date(sub("bigResultsDF_wEpitopes_(\\d{4}-\\d{2}-\\d{2})\\.csv", "\\1", basename(resultDFs))))]))
bigResults <- unique(bigResults[!is.na(bigResults$subtype),])

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
      runName <- paste0(subtype,"_epitopeAnalysis_t", threshold, "_", paste0(mySeason, collapse="."))
      print(runName)
      
      highFreqSites <- bigResults$highfreqMutSites[bigResults$subtype==subtype & 
                                                     bigResults$threshold==threshold & bigResults$mySeasons==mySeason]
      if(length(highFreqSites)>0){
        highFreqSites <- as.numeric(strsplit(highFreqSites, " ")[[1]])
        mutDF <- create_mut_df_fromList(bigMutDF, listMut = highFreqSites) %>% dplyr::rename(state=location) # subset every sequence to the set of hgih freq sites
        haplos <- mutDF[,c("label","date", "season", "UMich", "HAH3muts", "state", "highfreqMut" )]
        write.csv(haplos, paste0("haplos_perSeason/",runName,"_haplos.csv"), row.names = F)
      }
    }
  }
}




###############
# 3. Haplotype stats - from script 3
##############



rootDir <- paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/epitopeAnalysis/")
setwd(rootDir)

mySeasons <- unique(bigResults$mySeasons)


haploStatsDF <- data.frame(subtype =NA, threshold=0, season = NA, 
                           p90noMich=0,p50noMich=0,topXinMich_noMich = 0,ofTop50_nMissingInMich_noMich = 0,ofTop90_nMissingInMich_noMich = 0,
                           sumPropStrains_inMich_noMich = 0,sumPropStrains_missingMich_noMich = 0,
                           p90Mich=0, p50Mich=0, topXinMich_wMich = 0, ofTop50_nMissingInMich_wMich = 0, 
                           ofTop90_nMissingInMich_wMich = 0, sumPropStrains_inMich_wMich = 0, sumPropStrains_missingMich_wMich = 0)
i <- 0
for (subtype in subtypes){
  
  setwd(rootDir)
  
  for (threshold in thresholds){
    for (mySeason in mySeasons){
      
      i <- i+1
      runName <- paste0(subtype,"_epitopeAnalysis_t", threshold, "_", paste0(mySeason, collapse="."))
      print(runName)
      
      haploDF_noMich <- read.csv(paste0("haploDF_noMich_", runName,".csv"))
      haploDF_wMich <- read.csv(paste0("haploDF_wMich_", runName,".csv"))
      
      #update with haplos that were detected at a later date
      # identify the correct file
     # if(mySeason != "2021.2022.2023.2024"){
        
        fullHaploinfo_fileName <- paste0("haplos_perSeason/",runName , "_haplos.csv")
        fullHaploinfo <- read.csv(fullHaploinfo_fileName)
        
        haploDF_wMich2 <- fill_in_earlierHaplos(haploDF_wMich, fullHaploinfo, mySeason=mySeason)
        haploDF_noMich2 <- fill_in_earlierHaplos(haploDF_noMich, fullHaploinfo, mySeason=mySeason)

       
     # } else {
      #  haploDF_wMich2 <- haploDF_wMich
       # haploDF_noMich2 <- haploDF_noMich
      #}
      
      write.csv(haploDF_wMich2, paste0(rootDir,  "/haploDF_wMich_", runName,"_v2.csv"), row.names = F)
      write.csv(haploDF_noMich2, paste0(rootDir, "/haploDF_noMich_", runName,"_v2.csv"), row.names = F)
      
      # what would be more informative, below, than the number of syrains not in michigan, is the proportion that they represent
      # X% of strains have a Michigan haplo match, X% f srains do not
      
      haploStatsDF[i,] <-  c(subtype, threshold, paste(mySeason, collapse = "."),
                             unlist(haplotype_stats(haploDF_noMich2)), unlist(haplotype_stats(haploDF_wMich2)))
    }
  }
}

write.csv(haploStatsDF, paste0(rootDir, "bigHaploStatsDF_wEpitopes_", Sys.Date(), ".csv"))



############
# change in delays and prop detected
###########

#want to show X axis : proportion sampled
#y axis : number of NAs, delay in detection



rootDir <- "/Users/manonragonnet/Library/CloudStorage/OneDrive-Personal/Projects/LauringLab/Representativeness/haplo_stats3/epitopeAnalysis/"
setwd(rootDir)

myFiles <- list.files(pattern = "haploDF_noMich.*_v2", recursive = TRUE)

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
dropLines1 <- which(all_data$subtype=="H1" & all_data$season=="2021")
dropLines2 <- which(is.na(all_data$firstUSDate))  # weird result here with H1 2021 0.005 to check # fixed, should be empty

all_data <- all_data[- c(dropLines1, dropLines2),]
all_data$WEEK <- MMWRweek(all_data$firstUSDate)$MMWRweek
table(all_data$firstUS_season, all_data$subtype, all_data$threshold)

write.csv(all_data, paste0(rootDir, "/haploStats_wEpitopes_", Sys.Date(), ".csv"), row.names = F)






