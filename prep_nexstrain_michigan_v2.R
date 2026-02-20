# ============================================================================
# prep_nexstrain_michigan_v2.R
# Data preparation: combine GISAID and GenBank sequences. Deduplicates
# entries, assigns US state labels and Michigan/non-Michigan flags, and
# produces Nextstrain input files.
#
# Inputs: *_HA_genbank.csv, *_HA_gisaid.csv, *_HA_genbank_gisaid.fas
# Outputs: *_nextstrain.epi.tsv, *_nextstrain.fas
# Dependencies: cleanUS_states.R
# ============================================================================

library(ape)
library(dplyr)
library(lubridate)
library(gtools)


source("/Users/manonragonnet/OneDrive/Github/Michigan_rep/R_Michigan/cleanUS_states.R")
source("/Users/manonragonnet/OneDrive/Github/Michigan_rep/R_Michigan/functions_pick_strains_functions.R")

setwd("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/raw_data_Oct1/") 

tidyDates <- function(csv, dateCol="date"){ # also drops rows with incomplete dates
  dateColID <- which(names(csv)==dateCol)   
  csv$nextstrain_date <- as.Date(parse_date_time(csv[,dateColID], orders = "ymd"))
  csv <- csv%>% filter(!is.na(nextstrain_date))
  
  return(csv)
}

combineCSVs <- function(subtype){
  
  genbank_csv<- read.csv(paste0(subtype,"_HA_genbank.csv"))
  gisaid_csv <- read.csv(paste0(subtype,"_HA_gisaid.csv"))
  
  # somr print checks
  print(table(gisaid_csv$Subtype))
  
  gisaid_csv2 <- gisaid_csv%>%select(Isolate_Id, Isolate_Name, Collection_Date, Clade, Location, Host)# HA.Segment_Id, Originating_Lab,Genotype,
  genbank_csv2 <- genbank_csv%>%select(accession, name, date,  region, host) %>% #lineage,
    dplyr::rename(Location = region, Isolate_Id = accession, Host=host) 
  
  #tidy Dates
  genbank_csv2<- tidyDates(genbank_csv2, dateCol = "date")
  gisaid_csv2<- tidyDates(gisaid_csv2, dateCol = "Collection_Date")
  
  #add seqNames
  gisaid_csv2 <-  gisaid_csv2 %>% mutate(seqName = paste0(Isolate_Id, "|", Isolate_Name, "|", nextstrain_date))
  genbank_csv2 <- genbank_csv2 %>% mutate(Isolate_Name = sub(".*\\((.*)\\).*", "\\1", name)) %>%
    mutate(seqName = paste0(Isolate_Id, "|", Isolate_Name, "|", nextstrain_date))
  
  # relabel columns
  genbank_csv2 <- genbank_csv2%>% select(-c(name, date)) %>% mutate(Clade=NA)%>% mutate(seqSource = "genbank")
  gisaid_csv2 <- gisaid_csv2%>% select(-c( Collection_Date))%>% mutate(seqSource = "gisaid")
  
  # concat files
  csv <- smartbind(genbank_csv2, gisaid_csv2) 
  
  #1 assign all states
  #2 michigan/ non-Michigan
  csv <- modifyMich(csv) %>% select(-c(state1, state2))
  
  cat("[prep] ", subtype, " | after modifyMich | total:", nrow(csv),
      "| UMich:", sum(csv$UMich == "UMich"), "| notUOM:", sum(csv$UMich == "notUOM"), "\n")
  csv$subtype <- subtype
  csv$seqName <- gsub(" | ", "|", csv$seqName, fixed = T)
  csv$seqName <- gsub("  |", "|", csv$seqName, fixed = T)
  csv$seqName <- gsub(" |", "|", csv$seqName, fixed = T)
  csv$seqName <- gsub("| ", "|", csv$seqName, fixed = T)
  csv$seqName <- gsub(" ", "_", csv$seqName)
  
 
  # print out numbers by state
  
  print(table(csv$state))
  
  #3 deduplication round 1 - only keep 1 per isolate name
  # my later deduplication script relies on aligned sequences, so it has to be done after nextstrain
  
  csv2<- csv %>%
    # Give priority score: gisaid = 1, genbank = 2
    mutate(source_rank = if_else(seqSource == "gisaid", 1L, 2L)) %>%
    arrange(Isolate_Name, source_rank) %>%   # ensure gisaid comes first
    distinct(Isolate_Name, .keep_all = TRUE) %>%  # keep first per isolate
    select(-source_rank)
  
  print(nrow(csv))
  print(paste("deduped:", nrow(csv2)))
  cat("[prep] ", subtype, " | after dedup | total:", nrow(csv2),
      "| UMich:", sum(csv2$UMich == "UMich"), "| notUOM:", sum(csv2$UMich == "notUOM"), "\n")
  
  write.csv(csv2,paste0(subtype, "_GISAID_GENBANK_unique.csv" ), row.names = F)
  
}


create_nextstrain_epi <- function(subtype){
  #lat_long <- read.csv(paste0(topfolder, "allCountries_lat_longs.csv"))
  csv <- read.csv(paste0(subtype, "_GISAID_GENBANK_unique.csv" ))
  fasta <- read.FASTA(paste0(paste0( subtype, "_HA_genbank_gisaid.fas")))
  names(fasta) <- gsub(" ", "_", names(fasta))
  seqNames <- names(fasta)
  
  csv <-  csv %>%    mutate(strain = seqName)
  
  csv_nextstrain <- csv
  csv_nextstrain$date <- csv_nextstrain$nextstrain_date
  csv_nextstrain$virus <- "influenza"
  csv_nextstrain$location <- csv_nextstrain$state
  
  # keep overlap between epi file and fasta file only
  nrow(csv_nextstrain)
  matches <- intersect(csv_nextstrain$strain, seqNames)
  csv_nextstrain$strain[which(!csv_nextstrain$strain %in% seqNames)]
  
  epi <- csv_nextstrain[csv_nextstrain$strain %in% seqNames,]
  fasta2 <- fasta[seqNames %in% csv_nextstrain$strain]
  fasta2 <- fasta2[unique(names(fasta2))]
  print(nrow(csv_nextstrain))
  print(nrow(epi))
  print(length(fasta2))
  
  epi <- epi[,c("strain", "date", "virus",  "location",  "UMich")] #"Continent",
  
  #return(epi, fasta2)
  
  write.table(epi, file=paste0("HA.",subtype, "_nextstrain.epi.tsv"), quote=FALSE, sep='\t', row.names = F)
  write.FASTA(fasta2, file=paste0("HA.",subtype, "_nextstrain.fas"))
  
}






combineCSVs("H1")
combineCSVs("H3")

create_nextstrain_epi("H1")
create_nextstrain_epi("H3")
