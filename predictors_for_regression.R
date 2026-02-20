# ============================================================================
# predictors_for_regression.R
# Prepare external predictor variables for regression. Processes ILINet
# surveillance data, US Census population estimates, geographic distances
# to Michigan, and per-state sequencing counts.
#
# Inputs: ILINet.csv, *_ALLmuts.rds
# Outputs: episize_stats.csv, state_metrics.csv, state_seq_numbers.csv
# Dependencies: functions_pick_strains_functions.R
# ============================================================================

# create the  big regression dataframe

library(dplyr)
library(tidycensus)
library(sf)
library(tigris)



setwd("/Users/manonragonnet/Library/CloudStorage/OneDrive-Personal/Projects/LauringLab/Representativeness/")
source("/Users/manonragonnet/Library/CloudStorage/OneDrive-Personal/Github/Michigan_rep/R_Michigan/Representativeness/functions_pick_strains_functions.R")

#what I want to add to the dataframe
# I need to generate some daily epi counts of epidemics size

ILInet <- read.csv("raw-data/ILINet.csv", skip=1)

ILInet <- ILInet %>% select(YEAR, WEEK, ILITOTAL) %>%
  mutate(SEASON = ifelse(WEEK>=40, YEAR, YEAR-1))

seasonTotals <- ILInet %>% dplyr::group_by(SEASON) %>%
  dplyr::summarise(SEASON_TOTAL = sum(ILITOTAL, na.rm = TRUE)) 

ILInet <- merge(ILInet, seasonTotals, by="SEASON", all=T) %>%
  dplyr::arrange(YEAR, WEEK) %>%
  dplyr::group_by(SEASON) %>%
  dplyr::mutate(cumILI = cumsum(ILITOTAL)) %>%
  dplyr::mutate(EPIDEMIC_PROP = cumILI/SEASON_TOTAL)%>% ##--> I will use season total and epidemic prop in my model
  dplyr::ungroup() %>%
  dplyr::mutate(FLUWEEK = ifelse(WEEK>=40, WEEK-39, WEEK+14))


write.csv(ILInet, "episize_stats.csv", row.names = F)

# state sizes
# install.packages("tidycensus")

# Set your Census API key - olny once
#census_api_key("YOUR_API_KEY", install = TRUE)

# Get most recent population estimates
pop_df <- get_estimates(geography = "state", product = "population", year = 2021) %>%
  filter(variable=="POPESTIMATE") %>%  select(NAME, value)


# distance to Mchigan
# Download simplified US states polygons
states <- states(cb = TRUE)

# Calculate state centroids
centroids <- st_centroid(states)

# Get distance matrix between state centroids (in meters)
dist_matrix <- st_distance(centroids)

# Convert to kilometers
dist_km <- units::set_units(dist_matrix, "km")

# Label rows/columns
state_names <- states$NAME
dimnames(dist_km) <- list(state_names, state_names)
dist_km <- as.data.frame(dist_km)
dist_km$state <- state_names
dist_to_Mich <- dist_km[, c("state","Michigan")]


state_metrics <- merge(pop_df, dist_to_Mich, by.x="NAME", by.y = "state")
state_metrics$Michigan <- as.numeric(state_metrics$Michigan)
names(state_metrics) <- c("state", "pop_size", "dist_to_Mich")

write.csv(state_metrics, "state_metrics.csv", row.names = F)


##############
# seasonal variability

#subtypes <- c("H1", "H3")
#mySeasons <- c(2021, 2022, 2023, 2024)

#df <- data.frame(subtype = c(rep("H1",4), rep("H3", 4)), season = rep(mySeasons,2), seasonal_entropy = NA)

#for (i in 1:nrow(df)){
 # print(i)
#  df$seasonal_entropy[i] <- seasonal_entropy(df$subtype[i], df$season[i])
#}




# number fo sequences per state/ season
H1 <- readRDS(paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/H1_ALLmuts.rds"))
H3 <- readRDS(paste0("/Users/manonragonnet/OneDrive/Projects/LauringLab/Representativeness/haplo_stats3/H3_ALLmuts.rds"))

H1 <- H1 %>% select(c(label, location, season)) %>% mutate(subtype = "H1")
H3 <- H3 %>% select(c(label, location, season)) %>% mutate(subtype = "H3")

H1H3 <- rbind(H1, H3)
df <- as.data.frame(table(H1H3$location, H1H3$season))
names(df) <- c("firstState", "season", "nSeqSubmitted")

write.csv(df, "state_seq_numbers.csv", row.names = F)


