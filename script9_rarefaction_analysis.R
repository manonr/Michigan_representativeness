# ============================================================================
# rarefaction_analysis.R
# Rarefaction/coverage analysis of haplotype diversity using iNEXT.
# Estimates sample coverage per state and the number of sequences needed
# for 95% coverage.
#
# Inputs: haploDF_*.csv
# Outputs: SupFig2 rarefaction figure, TableS2/S3
# Dependencies: state_analysis_functions.R, iNEXT
#
# KNOWN ISSUES:
#   - Lines ~194-232: dead exploratory code (iNEXT tests on specific files)
# ============================================================================

#Rarifcatio analysis of haplotypes


library(devtools) 
library(iNEXT) 
library(ggplot2)

source("/Users/manonragonnet/OneDrive/Github/Michigan_rep/R_Michigan/Representativeness/state_analysis_functions.R")

setwd("/Users/manonragonnet/Library/CloudStorage/OneDrive-Personal/Projects/LauringLab/Representativeness/haplo_stats3")


seasons <- c(2021:2024, "2021.2022.2023.2024")
subtypes <- c("H3", "H1")
thresholds <- c(0.1)#, 0.01,0.005, 0.05, 0.15,0.2)
states <- unfactor(stateCounts_df$state)

df <- data.frame(subtype = "H3", state = "Michigan", season = "2021", threshold=0.1, nSeq  = 0, estimateD = 0 )
df2 <- data.frame(subtype = "H3",  season = "2021", threshold=0.1, nSeq  = 0, estimateD = 0 )
dfline <- 0 
df2line <- 0 
for (subtype in subtypes){
  for (threshold in thresholds){
    for (season in seasons){
      df2line <- df2line+1
      full_distr<- read.csv(paste0("haploDF_wMich_", subtype,"_", threshold,"_",  season, ".csv"))
      
      est_95 <- estimateD(
        full_distr$Count,
        datatype = "abundance",
        base = "coverage",
        level = 0.95
      )  
      
      cat("[script9] ", subtype, " | season", season, "| threshold", threshold,
          "| national haplos:", nrow(full_distr), "| national seqs:", sum(full_distr$Count), "\n")
      df2[df2line,] <- c(subtype, season,threshold, sum(full_distr$Count), round(est_95$m[1]))
      
      #setwd("stateProxyAnalysis/")
      for (state in states){
        print(paste(subtype, state, season, threshold))
        fileName <- paste0("stateProxyAnalysis/haploDF_", subtype, "_", state, "_stateProxyAnalysis_t", threshold, "_", season, "_v2.csv")
        if (file.exists(fileName)){
          dfline <- dfline+1
          state_distr <- read.csv(fileName)
          
          full_distr <- full_distr[,1:2]
          names(full_distr)[2] <- paste0(names(full_distr)[2], "_US")
          state_distr <- state_distr[,1:2]
          names(state_distr)[2] <- paste0(names(state_distr)[2], "_NoState")
          haplo_distr <- merge(full_distr, state_distr, by="Haplotype")
          haplo_distr$count_state <- haplo_distr$Count_US - haplo_distr$Count_NoState
          cat("[script9] ", subtype, " |", state, "| season", season,
              "| state seqs:", sum(haplo_distr$count_state),
              "| national seqs:", sum(haplo_distr$Count_US), "\n")
          #haplo_distr$Count_NoState <-NULL
          
          # out <- iNEXT(
          #    x = list(
          #      State = haplo_distr$count_state,
          #     National = haplo_distr$Count_US
          #   ),
          #    q = 0,
          #    datatype = "abundance"
          #  )
          
          #ggiNEXT(out, type = 2)
          est <- estimateD(haplo_distr$Count_US, datatype = "abundance", base = "size", level =sum(haplo_distr$count_state))
          
          df[dfline, ]<- c(subtype, state, season, threshold, sum(haplo_distr$count_state),est$SC[1])
          
        }
      }
    }
  }
}

df <- df[order(df$estimateD),]
df2 <- df2[order(df2$estimateD),]
df$season <- as.character(df$season)
df2$season <- as.character(df2$season)
df$subtype <- as.character(df$subtype)
df2$subtype <- as.character(df2$subtype)
df$estimateD <- as.numeric(df$estimateD)
df2$estimateD <- as.numeric(df2$estimateD)



# plot

library(dplyr)
library(ggplot2)
library(ggrepel)


# (optional but recommended) make sure types are right
df_cov <- df %>%
  mutate(
    nSeq = as.numeric(nSeq),
    estimateD = as.numeric(estimateD),
    season = as.character(season),
    subtype = as.character(subtype),
    state = as.character(state)
  )



# If your seasons are numeric like 2021, 2022, etc., keep only those
df1 <- df_cov %>%
  filter(season %in% c("2021", "2022", "2023", "2024")) %>%
  mutate(
    season = factor(season, levels = c("2021", "2022", "2023", "2024")),
    subtype = factor(subtype)  # or set levels = c("H1","H3") if you want
  )

df_labels <- df1 %>%
  group_by(subtype, season) %>%
  filter(nSeq > 0) %>%   # optional: ignore zero-sample states
  arrange(estimateD) %>%
  slice_head(n = 3) %>%  # bottom 3
  bind_rows(
    df1 %>%
      group_by(subtype, season) %>%
      filter(nSeq > 0) %>%
      arrange(desc(estimateD)) %>%
      slice_head(n = 3)         # top 3
  ) %>%
  ungroup()

p1 <- ggplot(df1, aes(x = nSeq, y = estimateD)) +
  geom_point(size = 2, alpha = 0.8) +
  
  # RED vertical line: sequences needed for 95% SC
  geom_vline(
    data = df2%>% filter(!season %in% c("2021.2022.2023.2024")),
    aes(xintercept = estimateD),
    colour = "firebrick",
    linewidth = 1,
    alpha = 0.8
  ) +
  
  geom_text_repel(
    data = df_labels,
    aes(label = state),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0
  ) +
  
  facet_grid(rows = vars(season), cols = vars(subtype), scales = "free_x") +
  
  labs(
    x = "Number of sequences",
    y = "Estimated coverage of national diversity (per state)"
  ) +
  
  theme_minimal(base_size = 12)+ 
  scale_y_continuous(limits = c(0, 1), labels = scales::percent)

p1

ggsave("../figures/SupFig2_Rarefaction_fig_states.jpg",p1, height=6, width=7)




df_wide <- df1 %>%
  mutate(
    subtype_season = paste0(subtype, "_", season),
    estimateD = round(estimateD,2)
  ) %>%
  select(state, subtype_season, estimateD) %>%
  pivot_wider(
    names_from  = subtype_season,
    values_from = estimateD
  ) %>%
  mutate(across(-state, ~ as.numeric(as.character(.x)))) %>%
  mutate(
    mean_estimateD = round(rowMeans(
      select(., -state, -H1_2021),
      na.rm = TRUE
    ),2)
  ) %>%
  arrange(desc(mean_estimateD))


df_wide

write.csv(df_wide, "../figures/TableS2_estimateD_perSeasonSubtypeState2.csv", row.names = F)

write.csv(df2, "../figures/TableS3_nSeq_for95pGD.csv", row.names = F)





# not needed?


df <- read.csv("haploDF_noMich_H1_0.01_2024_v2.csv")
df

df2 <- df[, 2]

mb <- iNEXT(df2, q=0, datatype="incidence_freq", size=1000, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=50)

ggiNEXT(mb, type=1, se=TRUE, facet.var="None", color.var="Assemblage", grey=FALSE)  


out <- iNEXT(df2, q=c(0, 1, 2), datatype="abundance", endpoint=sum(df2)) 
#ggiNEXT(out, type=1, facet.var="Assemblage")

ggiNEXT(out, type=2, facet.var="None", color.var="Assemblage")



#######

df <- read.csv("test2.csv")

freq_state <- df[,3]
freq_nat <- df[,2]

out <- iNEXT(
  x = list(
    State = freq_state,
    National = freq_nat
  ),
  q = 0,
  datatype = "abundance"
)

ggiNEXT(out, type = 2)
estimateD(freq_nat, datatype = "abundance", base = "size", level =sum(freq_state))

