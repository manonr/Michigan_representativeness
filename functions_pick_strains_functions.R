# ============================================================================
# functions_pick_strains_functions.R
# Core utility functions for the haplotype representativeness analysis.
#
# Provides functions to:
#   - Read and translate influenza HA sequences (DNA/AA FASTA)
#   - Extract amino acid haplotypes at high-frequency mutation sites
#   - Build ranked haplotype frequency tables
#   - Track first detection dates (US-wide and Michigan/state-specific)
#   - Compute Shannon entropy of haplotype distributions
#   - Summarize detection delays and haplotype coverage statistics
#   - Downsample Michigan sequences for sensitivity analysis
#
# Dependencies: cleanUS_states.R (for extract_us_state used in modifyMich)
#
# M Ragonnet-Cronin, 2024
# ============================================================================

library(lubridate)
library(tidyr)
library(dplyr)
require(ggplot2)
library(bioseq)
library(ape)
library(stringr)
library(ggpubr)
library(readr)
library(reshape2)
library(tibble)
library(Biostrings)


source("/Users/manonragonnet/OneDrive/Github/Michigan_rep/R_Michigan/cleanUS_states.R")

transform_auspice_df <- function(auspice_file = "auspicedf.json.tsv"){
  
  auspice_df <- read.table(auspice_file, sep="\t", header = T)
  auspice_df$HAmut <- unlist(lapply(auspice_df$mutations, function(x) {strsplit(x, "HA:")[[1]][2]}))
  auspice_df$HAmut <- gsub("\\[|\\]| |\\{|\\}", "", auspice_df$HAmut)
  #
  auspice_df$year <- year(unlist(lapply(auspice_df$num_date,date_decimal)))
  auspice_df$month<- month(as.Date(do.call("c",(lapply(lapply(auspice_df$num_date,date_decimal), function(x) {as.Date(x)})))))
  for (i in 1:nrow(auspice_df)){
    if(auspice_df$month[i]<10){
      auspice_df$season[i] <-as.numeric(auspice_df$year[i])-1
    } else {
      auspice_df$season[i] <-as.numeric(auspice_df$year[i])
    }
  }
  names(auspice_df)[names(auspice_df)=="name"] <- "label"; auspice_df <- as_tibble(auspice_df)
  return(auspice_df)
}

auspice_getMut <- function(auspice_df){
  HAmut <- gsub("NA,|,NA", "",paste(auspice_df$HAmut, collapse=","))
  HAmut <- unique(strsplit(HAmut, ","))
  HAmutSites <- sort(as.numeric(unique(unlist(lapply(HAmut, parse_number)))))
  #HAmut[as.numeric(unlist(order(HAmutSites)))]
  return(list(HAmut, HAmutSites))
}


make_mut_df <- function(dna_fasta, 
                        listMut = c(3,4,5,8,9,10,11,12,13,14,15,16,18,19,21,22,23,24,25,26,28,31,34,35,36,37,40,41,43,46,47,48,49,50,53,59,60,61,62,64,65,66,69,70,72,73,74,76,78,79,81,82,83,85,89,91,94,95,96,97,98,99,101,104,107,108,110,112,117,119,120,121,122,125,128,130,135,137,138,140,141,142,144,145,147,149,151,153,154,156,158,159,160,161,162,166,171,172,173,174,175,176,177,179,180,182,183,184,185,187,188,189,190,193,195,198,199,202,204,205,206,208,209,210,212,213,214,215,217,218,219,220,223,224,225,228,229,230,233,235,238,239,241,242,243,245,246,248,249,258,260,264,275,276,277,278,280,283,285,287,288,290,291,292,294,296,299,305,307,313,315,316,320,323,325,326,327,328,339,342,344,347,349,362,363,377,388,391,395,396,397,400,402,406,418,422,424,427,429,431,434,435,443,446,451,454,456,457,458,466,468,469,478,484,488,490,492,494,495,496,500,503,505,506,509,511,517,519,520,521,524,529,531,532,538,540,541,545,546,547,552,557,563,566)){
  #aa_seq <- bioseq::read_fasta("aligned_aa_Michigan_plus_vax_H3N2_HA.fasta", type="AA")
  nt_seq <- bioseq::read_fasta(dna_fasta, type="DNA")
  #vax_strains <- c("A/Singapore/INFIMH-16-0019/2016")
  
  
  names(nt_seq) <- gsub(" ", "", gsub("\r", "", names(nt_seq)))
  fra_data <- tibble(label = names( nt_seq ), sequence =  nt_seq )
  fra_data$category <- "US"
  fra_data$category[which(grepl("Michigan", fra_data$label))] <- "Michigan"
  fra_data$category[grep("NODE_", fra_data$label)] <- "node_seq"
  #fra_data$category[which(!grepl("NODE_|Michigan", fra_data$label))] <- "vax_seq"
  #fra_data$label[match(vax_seq$V1, fra_data$label)] <- "vax_seq"
  cds_start <- gregexpr(pattern ='ATG',fra_data$sequence[[1]])[[1]][1]
  cds_end <- cds_start +1698
  
  fra_data$cds <- unlist(lapply(fra_data$sequence, function(y) {bioseq::as_dna(substr(y, cds_start,cds_end ))}))
  fra_data$aa <- unlist(lapply(fra_data$cds,function(x) {bioseq::seq_translate(bioseq::as_dna(x), code = 1, codon_frame = 1, codon_init = T)}))
  
  # jesse person to person variation in targeted sites https://elifesciences.org/articles/49324
  # clade detrmining mutation in nextstrain https://github.com/nextstrain/seasonal-flu/blob/master/config/h3n2/ha/clades-long.tsv
  
  for(i in listMut){
    fra_data$newCol <- unlist(lapply(fra_data$aa, function(y) { paste(lapply( i, function(x) {substr(y, x,x )}), collapse="")}))
    names(fra_data)[names(fra_data)=="newCol"] <- paste0("aapos", i)
  }
  
  fra_data$HAH3muts <- unlist(lapply(fra_data$aa, function(y) { paste(lapply(listMut, function(x) {substr(y, x,x )}), collapse="")}))
  
  
  return(fra_data)
  
  
}


fill_in_earlierHaplos<- function(df, fullHaploinfo, mySeason){
  df$last_haplo <- NA
  if (mySeason != "2021.2022.2023.2024"){
    df$USseason_dominant <- df$firstUS_season
  } else {
    df$USseason_dominant <- mySeason
  }
  for (i in 1:nrow(df)){
    #print(i)
    haploMatches <- fullHaploinfo%>% filter( HAH3muts==df$Haplotype[i])   # now for all rows, we replace first haplo with the first haplo EVER, rather than in that seaosn only
    firstMatch <- haploMatches[which(haploMatches$date==min(haploMatches$date)),][1,]
    
    df$firstUS_season[i] <- season_from_single_date(firstMatch$date)
    df$firstUSDate[i] <-firstMatch$date
    df$firstUSstrain[i] <- firstMatch$label
    df$firstState[i] <- firstMatch$state
    
    michMatches <- fullHaploinfo%>% filter( HAH3muts==df$Haplotype[i] & UMich=="UMich") # and even if a mich match was found, we replace it with the first mich haplo EVER
    
    if(nrow(michMatches)>0){
      firstMichMatch <- michMatches[which(michMatches$date==min(michMatches$date)),][1,]
      df$firstMichSeason[i] <- season_from_single_date(firstMichMatch$date)
      df$firstMichDate[i] <-firstMichMatch$date
      df$firstMichStrain[i] <- firstMichMatch$label
    }     else {
      df$firstMichSeason[i] <- NA;      df$firstMichDate[i] <-NA;      df$firstMichStrain[i] <- NA
    }
    df$last_haplo[i] <- max(fullHaploinfo$date[fullHaploinfo$HAH3muts==df$Haplotype[i]])
  }
  df$delay <- as.Date(df$firstMichDate)-as.Date(df$firstUSDate)
  df$delay[is.na(df$firstMichStrain)]<- NA
  return(df)
}

fill_in_earlierHaplos_allStates<- function(df, fullHaploinfo, mySeason, mystate){
  df$last_haplo <- NA
  if (mySeason != "2021.2022.2023.2024"){
    df$USseason_dominant <- df$firstUS_season
  } else {
    df$USseason_dominant <- mySeason
  }
  for (i in 1:nrow(df)){
    #print(i)
    haploMatches <- fullHaploinfo%>% filter( HAH3muts==df$Haplotype[i])   # now for all rows, we replace first haplo with the first haplo EVER, rather than in that seaosn only
    firstMatch <- haploMatches[which(haploMatches$date==min(haploMatches$date)),][1,]
    
    df$firstUS_season[i] <- season_from_single_date(firstMatch$date)
    df$firstUSDate[i] <-firstMatch$date
    df$firstUSstrain[i] <- firstMatch$label
    df$firstState[i] <- firstMatch$state
    
    stateMatches <- fullHaploinfo%>% dplyr::filter( HAH3muts==df$Haplotype[i] & state==mystate) # and even if a mich match was found, we replace it with the first mich haplo EVER
    
    if(nrow(stateMatches)>0){
      firstStateMatch <- stateMatches[which(stateMatches$date==min(stateMatches$date)),][1,]
      df$firstStateSeason[i] <- season_from_single_date(firstStateMatch$date)
      df$firstStateDate[i] <-firstStateMatch$date
      df$firstStateStrain[i] <- firstStateMatch$label
    } else {
      df$firstStateSeason[i] <- NA;      df$firstStateDate[i] <-NA;      df$firstStateStrain[i] <- NA
    }
    df$last_haplo[i] <- max(fullHaploinfo$date[fullHaploinfo$HAH3muts==df$Haplotype[i]])
  }
  df$delay <- as.Date(df$firstStateDate)-as.Date(df$firstUSDate)
  df$delay[is.na(df$firstStateStrain)]<- NA
  return(df)
}

season_from_single_date <- function(date1){
  if(!is.na(date1)){
    #print(date1)
    myYear <- year(date1)
    myMonth<- month(date1)
    if(myMonth<10){
      mySeason <-as.numeric(myYear)-1
    } else {
      mySeason <-as.numeric(myYear)
    }
    return(mySeason)
  } else {
    return(NA)}
}

getMutList <- function(aa_df, colname = "aapos", myCutoff=0.05){
  # need to put in an exception, if I have fewer than 100 sequences, .. use a hgiher threshold?
  print(unique(aa_df$season))
  
  aapos_cols <- grep("aapos", names(aa_df), value = TRUE)
  cols_with_na <- aapos_cols[colSums(is.na(aa_df[aapos_cols])) > 0]
  aa_df <- aa_df[, !names(aa_df) %in% cols_with_na]
  
  if(nrow(aa_df)>=20){   # else take no sites
    
    unique_haplo1<- aa_df[!duplicated(aa_df$aa),] # these are the unique sequences
    
    #if(nrow(unique_haplo1)<50){
    # myCutoff <- 0.3
    #  print("fewer than 50 unique haplotypes")
    #} else {
    #  myCutoff <- 0.1
    #  print(" 50 or more unique haplotypes")
    #}
    
    highfreqMutSites <- c()
    for (myCol in names(unique_haplo1)[grep(colname, names(unique_haplo1))]){
      #print(myCol)
      myaaprops <- table((unique_haplo1[,myCol]))
      otherProps <- sum(myaaprops[-unname(which(myaaprops==max(myaaprops)))])/sum(myaaprops)
      if (otherProps>myCutoff){
        highfreqMutSites <- c(highfreqMutSites, myCol)
        #print(myaaprops)
        #} else if(length(myaaprops)>3 &length(which(as.numeric(myaaprops)>0.05))==1){
        # print(myaaprops)
      }
      
    }
    highfreqMutSites <- as.numeric(gsub(colname, "", highfreqMutSites))
    return(highfreqMutSites)
  }  else {
    print("too few sequences to estimate high frequency mutation sites")
  }
}

getMutCounts <- function(aa_df, colname = "aapos", listSites=c(10,20)){
  
  aa_df_realseq <- aa_df[aa_df$category != "ancestor_seq" & aa_df$category != "node_seq",]
  mut_df_counts <- data.frame(site = listSites, nondom_freq = "unknown", domAA = "unknown", minAAs = "unknown")
  
  for (mySite in listSites){
    mySite2 <-  paste0(colname, mySite)
    print(mySite2)
    
    myCol <- match(mySite2, names(aa_df_realseq))
    
    if(is.na(myCol)){  # first check that we even have a column for that aa
      mut_df_counts$nondom_freq[ mut_df_counts$site==mySite] <- "no variation"
    } else {
      myaaprops <- table((aa_df_realseq[,mySite2]))
      otherProps <- sum(myaaprops[-unname(which(myaaprops==max(myaaprops)))])/sum(myaaprops)
      print(myaaprops)
      mut_df_counts$nondom_freq[ mut_df_counts$site==mySite] <- round(otherProps,2)
      mut_df_counts$domAA[ mut_df_counts$site==mySite] <- names(which(myaaprops==max(myaaprops)))
      mut_df_counts$minAAs[ mut_df_counts$site==mySite] <- paste(names(myaaprops[-unname(which(myaaprops==max(myaaprops)))]), collapse="-")
      print(paste("nondominant AA in",  round(otherProps*100,2), "% of sequences"))
    }
  }
  return( mut_df_counts)
} 

create_big_mut_df <- function(subtype, run){
  setwd(paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/Nextstrain/",run,  "/nextstrain_", subtype, "/"))
  auspice_df <- transform_auspice_df(paste0("HA_", subtype,"_nextstrain.json.tsv"))
  nextstrain_epi <- read.csv(paste0("data/HA.", subtype, "_nextstrain.epi.tsv"), sep="\t")
  #nextstrain_epi <- modifyMich(nextstrain_epi)
  
  
  # identify ALL mutations, and generate mutation frequencies
  HAmutSites <- auspice_getMut(auspice_df = auspice_df)[[2]]
  print(paste("all mutations: ",
              paste(HAmutSites, collapse=" ")))
  
  aa_df1 <- make_mut_df(paste0("results/aligned_nt_", subtype, "_nextstrain_%GENE.fasta"), listMut = HAmutSites) 
  aa_df1b <- merge( auspice_df[,c("label", "season")], aa_df1, by="label")
  
  USTable <- aa_df1b[aa_df1b$category !="ancestor_seq" & aa_df1b$category !="node_seq",]
  USTable <- merge(USTable, nextstrain_epi, by.x="label", by.y="strain", all.x=T)
  
  print(table(USTable$season)) 
  
  #save the big table
  return(USTable)
  
}


modifyMich <- function(nextstrain_epi){
  #modify as necessary
  nextstrain_epi <- nextstrain_epi
  
  #1 tidy up all state names
  nextstrain_epi$state1 <- extract_us_state(nextstrain_epi$Location)
  nextstrain_epi$state2 <- extract_us_state(nextstrain_epi$seqName)
  nextstrain_epi <- nextstrain_epi %>%
    mutate(state = coalesce(state1, state2))
  
  #1b extra step here want to change the state to Maryland for all those Baltimore sequences
  #print("TO DO!!! change the state to Maryland for all those Baltimore sequences")
  nextstrain_epi$state[grep("Baltimore", nextstrain_epi$seqName)] <- "Maryland"
  nextstrain_epi$state[grep("MICHGAN", nextstrain_epi$seqName)] <- "Michigan"
  nextstrain_epi$state[grep("Columbia", nextstrain_epi$seqName)] <- "Maryland"
  nextstrain_epi$state[grep("COLUMBIA", nextstrain_epi$seqName)] <- "Maryland"
  nextstrain_epi$state[grep("San Diego", nextstrain_epi$seqName)] <- "California"
  
  #2 find all UOM sequences and relabel UMich
  michList1 <- nextstrain_epi$seqName[grep("A/Michigan/CEIRS|A/Michigan/UM|A/Michigan/UOM|A/Michigan/UOMHS|A/Michigan/UOMMH|A/Michigan/ASJ|A/Michigan/UMHM|A/Michigan/MIS|A/Michigan/IVY|A/Michigan/HFH",
                                            nextstrain_epi$seqName)]  #
  #michList2 <- nextstrain_epi$strain[grep("Lauring", nextstrain_epi$strain)]
  
  
  nextstrain_epi$UMich[nextstrain_epi$seqName %in% michList1] <- "UMich"
  nextstrain_epi$UMich[!nextstrain_epi$seqName %in% michList1] <- "notUOM"
  
  #table(nextstrain_epi$UMich[grep("Michigan", nextstrain_epi$state)],unlist(lapply(nextstrain_epi$strain[grep("Michigan", nextstrain_epi$state)], function(x) {strsplit(strsplit(x, "A/Michigan/")[[1]][2], split = "[-0-9]", perl = TRUE)[[1]][1]})))
  
  
  return(nextstrain_epi)
  
} 

create_mut_df <- function(USTable, threshold, mySeason){
  
  
  print(unique(USTable$season)) # check that we are now only looking at the seasons of interest
  
  USTable_bySeason<- split.data.frame(USTable, USTable$season)
  
  # identify high freq mutations, in each season
  # I do think it makes sense to identify mutations using the Michigan sequences as well
  highfreqMutSites <- lapply(USTable_bySeason, function(x) {getMutList(x, "aapos", myCutoff=threshold)})
  if (length(unlist(highfreqMutSites))>0){
    ourYears <- which(as.numeric(names(highfreqMutSites))%in% mySeason)
    highfreqMutSites2 <- sort(unique(as.numeric(unlist(highfreqMutSites[ourYears]))[!is.na(as.numeric(unlist(highfreqMutSites[ourYears])))]))
    print(paste(" high freq mutations at threshold", threshold, ": ",
                paste(highfreqMutSites2, collapse=" ")))
    
    # make dataframe, for high freq mutations only
    # needs ot be changed to take big USdf
    # US table
    aaposCols <- paste0("aapos", highfreqMutSites2)
    non_aapos_cols <- colnames(USTable)[!grepl("aapos", colnames(USTable))]
    exact_matches <- colnames(USTable)[colnames(USTable) %in% c(aaposCols, non_aapos_cols)]
    USTable2 <- USTable[,exact_matches] %>% filter(season %in% mySeason)
    USTable2$HAH3muts <- apply(USTable2[, aaposCols, drop = FALSE], 1, paste0, collapse = "")
    USTable2$cds <- NULL; USTable2$virus <- NULL; USTable2$aa<- NULL; USTable2$sequence <- NULL
    USTable2$highfreqMut <- paste(highfreqMutSites2, collapse=" ")
    
    return(USTable2)
  } 
}

create_mut_df_fromList <- function(USTable, listMut){
  
  print(unique(USTable$season)) # check that we are now only looking at the seasons of interest
  print("no subsetting by season")
  
  # identify high freq mutations, in each season
  # I do think it makes sense to identify mutations using the Michigan sequences as well
  
  if (length(unlist(listMut))>0){
    
    # make dataframe, for high freq mutations only
    # needs ot be changed to take big USdf
    # US table
    aaposCols <- paste0("aapos", listMut)
    non_aapos_cols <- colnames(USTable)[!grepl("aapos", colnames(USTable))]
    exact_matches <- colnames(USTable)[colnames(USTable) %in% c(aaposCols, non_aapos_cols)]
    USTable2 <- USTable[,exact_matches] 
    USTable2$HAH3muts <- apply(USTable2[, aaposCols, drop = FALSE], 1, paste0, collapse = "")
    USTable2$cds <- NULL; USTable2$virus <- NULL; USTable2$aa<- NULL; USTable2$sequence <- NULL
    USTable2$highfreqMut <- paste(listMut, collapse=" ")
    
    return(USTable2)
  } 
}




make_haploDF <- function(haplos){
  haplos$cumsum<-cumsum(haplos$Freq) 
  haplos$propDivCov <- haplos$cumsum/sum(haplos$Freq) # this is the proportion of the diversity at these sites covered
  names(haplos) <- c("Haplotype", "Count", "Cumulative count", "Cumulative_proportion")
  haplos$frequency_order <- 1:nrow(haplos)
  return(haplos)
}

haploOrder <- function(USTable, inclMich=c(F, T, "only")){
  if(isFALSE(inclMich)){
    haplos <- as.data.frame(sort(table(USTable$HAH3muts[USTable$UMich=="notUOM"]), decreasing = T), stringsAsFactor=F)
    #drop haplos that have * or X
    badHaplos <- grep("\\*|X", haplos$Var1)
    if (length(badHaplos)>0) {haplos <- haplos[-badHaplos,]}
    haplos2 <- make_haploDF(haplos)
    
  } else if (isTRUE(inclMich)){
    haplos <- as.data.frame(sort(table(USTable$HAH3muts), decreasing = T), stringsAsFactor=F)
    badHaplos <- grep("\\*|X", haplos$Var1)
    if (length(badHaplos)>0) {haplos <- haplos[-badHaplos,]}
    haplos2 <- make_haploDF(haplos)
  } else if (inclMich=="only"){
    haplos <- as.data.frame(sort(table(USTable$HAH3muts[USTable$UMich=="UMich"]), decreasing = T), stringsAsFactor=F)
    #drop haplos that have * or X
    badHaplos <- grep("\\*|X", haplos$Var1)
    if (length(badHaplos)>0) {haplos <- haplos[-badHaplos,]}
    haplos2 <- make_haploDF(haplos)
  }
  
  michTable <- USTable[USTable$UMich=="UMich",]
  
  # sequences in Mich
  for (i in 1:nrow(haplos2)){
    # print(i)
    firstHaplo <- USTable[which(USTable$HAH3muts==haplos2$Haplotype[i] & 
                                  USTable$date==min(USTable$date[grep(haplos2$Haplotype[i], USTable$HAH3muts)])),][1,]
    #firstMichHaplo <- michTable[which(michTable$HAH3muts==haplos2$Haplotype[i] & 
    #                                   michTable$date==min(michTable$date[grep(haplos2$Haplotype[i], michTable$HAH3muts)])),][1,]
    
    firstMichHaplo <- michTable %>%
      dplyr::filter(HAH3muts == haplos2$Haplotype[i]) %>%
      dplyr::arrange(date) %>%
      dplyr::slice_head(n = 1)
    
    haplos2$firstUS_season[i] <- firstHaplo$season
    haplos2$firstUSDate[i] <-firstHaplo$date
    haplos2$firstState[i] <- firstHaplo$state
    haplos2$firstUSstrain[i] <- firstHaplo$label
    
    if(nrow(firstMichHaplo)>0){
      haplos2$firstMichSeason[i] <- firstMichHaplo$season
      haplos2$firstMichDate[i] <-firstMichHaplo$date
      haplos2$firstMichStrain[i] <- firstMichHaplo$label
    } else {
      haplos2$firstMichSeason[i] <- NA
      haplos2$firstMichDate[i] <- NA
      haplos2$firstMichStrain[i] <- NA
    }
    
    
    #haplos2$first_season[i] <- min(USTable$season[grep(haplos2$Haplotype[i], USTable$HAH3muts)])
    #haplos2$firstUSDate[i] <- min(USTable$date[grep(haplos2$Haplotype[i], USTable$HAH3muts)])
    #haplos2$firstMichDate[i] <- min(michTable$date[grep(haplos2$Haplotype[i], michTable$HAH3muts)])
  }
  haplos2$delay <- as.Date(haplos2$firstMichDate)-as.Date(haplos2$firstUSDate)
  haplos2$haplo_prop <- haplos2$Count/sum(haplos2$Count)
  
  return(haplos2)
  
}

haploOrder_pickState <- function(USTable, state = "Michigan", inclState=F){
  # first need to identify columns that are the state in question
  USTable$myState <- "no"
  USTable$myState[USTable$state==state] <- "yes"
  
  
  haplos <- as.data.frame(sort(table(USTable$HAH3muts[USTable$myState=="no"]), decreasing = T), stringsAsFactor=F)
  #drop haplos that have * or X
  badHaplos <- grep("\\*|X", haplos$Var1)
  if (length(badHaplos)>0) {haplos <- haplos[-badHaplos,]}
  haplos2 <- make_haploDF(haplos)
  
  stateTable <- USTable[USTable$myState=="yes",]
  
  for (i in 1:nrow(haplos2)){
    # print(i)
    firstHaplo <- USTable[which(USTable$HAH3muts==haplos2$Haplotype[i] & 
                                  USTable$date==min(USTable$date[grep(haplos2$Haplotype[i], USTable$HAH3muts)])),][1,]     
    
    firstStateHaplo <- stateTable %>%
      dplyr::filter(HAH3muts == haplos2$Haplotype[i]) %>%
      dplyr::arrange(date) %>%
      dplyr::slice_head(n = 1)
    
    haplos2$firstUS_season[i] <- firstHaplo$season
    haplos2$firstUSDate[i] <-firstHaplo$date
    haplos2$firstState[i] <- firstHaplo$state
    haplos2$firstUSstrain[i] <- firstHaplo$label
    
    if(nrow(firstStateHaplo)>0){
      haplos2$firstStateSeason[i] <- firstStateHaplo$season
      haplos2$firstStateDate[i] <-firstStateHaplo$date
      haplos2$firstStateStrain[i] <- firstStateHaplo$label
    } else {
      haplos2$firstStateSeason[i] <- NA
      haplos2$firstStateDate[i] <- NA
      haplos2$firstStateStrain[i] <- NA
    }
  }

  haplos2$delay <- as.Date(haplos2$firstStateDate)-as.Date(haplos2$firstUSDate)
  haplos2$haplo_prop <- haplos2$Count/sum(haplos2$Count)
  haplos2$thisState <- state
  
  return(haplos2)
  
}

haplotype_stats <- function(df){
  return(list(haplos90 = which(df$Cumulative_proportion>=0.9)[1],
              haplos50 = which(df$Cumulative_proportion>=0.5)[1],
              topXinMich = which(is.na(df$firstMichStrain))[1]-1,
              ofTop50_nMissingInMich= length(which(is.na(df$firstMichStrain[df$Cumulative_proportion<=0.5]))),
              ofTop90_nMissingInMich = length(which(is.na(df$firstMichStrain[df$Cumulative_proportion<=0.9]))),
              sumPropStrains_inMich= sum(df$haplo_prop[!is.na(df$firstMichStrain)]),
              sumPropStrains_missingMich = sum(df$haplo_prop[is.na(df$firstMichStrain)])))  
}
haplotype_stats_otherStates <- function(df){
  return(list(haplos90 = which(df$Cumulative_proportion>=0.9)[1],
              haplos50 = which(df$Cumulative_proportion>=0.5)[1],
              topXinState = which(is.na(df$firstStateStrain))[1]-1,
              ofTop50_nMissingInState= length(which(is.na(df$firstStateStrain[df$Cumulative_proportion<=0.5]))),
              ofTop90_nMissingInState = length(which(is.na(df$firstStateStrain[df$Cumulative_proportion<=0.9]))),
              sumPropStrains_inState= sum(df$haplo_prop[!is.na(df$firstStateStrain)]),
              sumPropStrains_missingState = sum(df$haplo_prop[is.na(df$firstStateStrain)])))  
}




# ---- helpers ----

shannon_entropy <- function(counts, base = 2) {
  counts <- counts[!is.na(counts)]
  counts <- counts[counts > 0]
  if (length(counts) == 0) return(NA_real_)
  
  p <- counts / sum(counts)
  if (base == 2) {
    -sum(p * log2(p))
  } else if (isTRUE(all.equal(base, exp(1)))) {
    -sum(p * log(p))
  } else {
    -sum(p * (log(p) / log(base)))
  }
}

summarise_entropy_by_group <- function(all_data,
                                       group_cols = c("season", "subtype", "threshold"),
                                       count_col = "Count",
                                       cum_col = "Cumulative.count",
                                       entropy_base = 2,
                                       tol = 0) {
  stopifnot(all(group_cols %in% names(all_data)))
  stopifnot(count_col %in% names(all_data))
  stopifnot(cum_col %in% names(all_data))
  
  # Requires dplyr
  dplyr::group_by(all_data, dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      n_haplotypes = sum(!is.na(.data[[count_col]])),
      sum_count = sum(.data[[count_col]], na.rm = TRUE),
      max_cumulative_count = suppressWarnings(max(.data[[cum_col]], na.rm = TRUE)),
      cumcount_equals_sumcount = isTRUE(
        abs(max_cumulative_count - sum_count) <= tol
      ),
      entropy = shannon_entropy(.data[[count_col]], base = entropy_base),
      effective_n = ifelse(is.na(entropy), NA_real_, entropy_base ^ entropy),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      max_cumulative_count = dplyr::if_else(is.infinite(max_cumulative_count), NA_real_, max_cumulative_count)
    )
}

# Optional: throw an error if any group fails the check

assert_cumcount_matches <- function(summary_df) {
  bad <- summary_df |>
    dplyr::filter(is.na(cumcount_equals_sumcount) | cumcount_equals_sumcount != TRUE)
  
  if (nrow(bad) > 0) {
    print(bad |> dplyr::slice_head(n = 10))
    stop(
      paste0(
        "Cumulative.count check failed or was NA for ", nrow(bad),
        " group(s). (Printed first 10 failing groups above.)"
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}






## NOTE: Old/unused functions have been moved to old_scripts/old_functions.R

# ---- Downsampling ----

drop_michProp <- function(mutDF, sample_prop){
  allUMich <- length(which(mutDF$UMich=="UMich"))
  sampleUMich_drop <- round((1-sample_prop)*allUMich,0)
  michToDrop <-sample(which(mutDF$UMich=="UMich"),sampleUMich_drop)
  if(length(michToDrop)>0){
    newmutDF <- mutDF[-michToDrop,]
  } else {
    newmutDF <- mutDF
  }
  return(newmutDF)
}

# ---- Detection delay summary ----

summarise_detection_delay <- function(df, top_prop = 0.5,
                                      detected_col = "detected",
                                      delay_col    = "delay",
                                      prop_col     = "haplo_prop") {

  df <- tibble::as_tibble(as.data.frame(df))

  df[[prop_col]]  <- as.numeric(df[[prop_col]])
  df[[delay_col]] <- as.numeric(df[[delay_col]])
  df[[detected_col]] <- as.integer(df[[detected_col]])

  summarise_one <- function(dat) {
    tibble::tibble(
      n_detected = sum(dat[[detected_col]] == 1L, na.rm = TRUE),
      n_not_detected = sum(dat[[detected_col]] == 0L, na.rm = TRUE),
      mean_delay_detected = mean(dat[[delay_col]][dat[[detected_col]] == 1L], na.rm = TRUE)
    )
  }

  df %>%
    dplyr::group_by(.data$subtype, .data$season, .data$threshold) %>%
    dplyr::group_modify(~{
      dat <- .x

      all_stats <- summarise_one(dat) %>% dplyr::mutate(set = "all")

      top_dat <- dat %>%
        dplyr::arrange(dplyr::desc(.data[[prop_col]])) %>%
        dplyr::mutate(cum_prop = cumsum(.data[[prop_col]])) %>%
        dplyr::filter(cum_prop <= top_prop)

      if (nrow(top_dat) == 0) {
        top_dat <- dat %>%
          dplyr::arrange(dplyr::desc(.data[[prop_col]])) %>%
          dplyr::slice_head(n = 1)
      }

      top_stats <- summarise_one(top_dat) %>%
        dplyr::mutate(set = paste0("top_", top_prop * 100, "pct"))

      dplyr::bind_rows(all_stats, top_stats)
    }) %>%
    dplyr::ungroup()
}

# ---- Scaling ----

scale_columns <- function(df, cols) {
  df %>%
    dplyr::mutate(
      dplyr::across(dplyr::all_of(cols), ~ as.numeric(scale(.x)))
    )
}
