# Name: Yu-Jing Chen
# Date: Jun 19th, 2025
# Task: Eviction Lab wk3 Practice 2
###GOAL: 
#Create an array that represents the Philadelphia MSA
#Suppose we are focusing on year 2022,which is the newest data available.
###

#install.packages("segregation")
#install.packages("remotes")
#remotes::install_github("arthurgailes/rsegregation")
library(rsegregation)
library(segregation)
library(tidyverse)
library(tidycensus)
library(here)


#1. Pull Data 
philly_msa <-tibble(
  state = c("PA", "PA", "PA", "PA", "PA", 
            "NJ", "NJ", "NJ", "NJ",
            "MD", "DE"),
  county = c("Bucks", "Chester", "Delaware", "Montgomery", "Philadelphia",
             "Burlington", "Camden", "Gloucester", "Salem",
             "Cecil", "New Castle")
)
variables <- c(
  total = "B25003_001",
  owner = "B25003_002",
  renter = "B25003_003"
)

get_data <-function(state, county) {
  get_acs(
    geography = "tract",
    state = state,
    county = county,
    variables = variables,
    year = 2022,
    survey = "acs5",
    output = "wide"
) %>% filter (
  totalE>0
) %>% select(
  GEOID, 
  totalE,
  ownerE,
  renterE
)
}

PM_data <- map2_dfr(
  philly_msa$state,
  philly_msa$county,
  get_data
)

#Save and read data locally
write.csv(PM_data, here("Philly_MSA_data.csv"),row.names = FALSE) 
PM_Data <- read.csv(here("Philly_MSA_data.csv"))

#2. Data Cleaning/Organizing (useful for divergence)
#Aim: Calculate proportions and total numbers for each group
PM_data <- PM_data %>% mutate(
  percent_renters = (renterE/totalE),
  percent_owners = (ownerE/totalE)
)
total_tracts <- nrow(PM_data)
#msa-level proportion
total_msa <- sum(PM_data$totalE)
renter_msa <- sum(PM_data$renterE)
owner_msa <- sum(PM_data$ownerE)

#3. Computing the 4 indices 
#3.1 Reshape data frame into long format (necessary for compuational needs)
long_pm <- pivot_longer(PM_data, cols = c(ownerE, renterE),
                        names_to = c("tenure"),
                        values_to = "n"
)
#3.2 Dissimilarity
D_tenure <-segregation::dissimilarity(
 long_pm,
 group = "tenure",
 unit = "GEOID",
 weight = "n"
)

#3.3 Isolation
I_tenure <-segregation::isolation(
  long_pm,
  group = "tenure",
  unit = "GEOID",
  weight = "n"
)

#3.4 Divergence Index
comparison_tenure <- c(
  renter = renter_msa / total_msa,
  owner = owner_msa/ total_msa
)
Diver_tenure <-divergence(
  PM_data$percent_renters,
  PM_data$percent_owners,
  population = PM_data$totalE,
  comparison = comparison_tenure,
  summed = TRUE
)

#3.5 MSA Entropy for Theil's H
entropy_msa <-segregation::entropy(
  long_pm,
  group = "tenure",
  weight = "n",
  base = 2
)
#for tract-level entropy
entropy_tract <- long_pm %>% 
  group_by(GEOID) %>%
  summarise(
    t_i = sum(n),
    E_i = segregation::entropy(pick(everything()),group = "tenure", weight ="n", base = 2)
  )

sigmpa_part <- sum(entropy_tract$t_i/total_msa*entropy_tract$E_i)

#plug all partd into Theil Index's formula
H_tenure <-(entropy_msa - sigmpa_part) / entropy_msa

#Defining a summary table
summary <-tibble(
  indices = c("Dissimilarity", "Isolation-renter", "Isolation-owner", "Divergence", "Theil's H"),
  MSA_values = c(0.408, 0.475, 0.740, 0.117, 0.184) 
  #rounded to 3 d.p.; un-comment the print statements below to view the full values. 
)
#print(D_tenure)
#print(I_tenure)
#print(Diver_tenure)
#print(H_tenure)

#4. Extension: macro-level segregation
#Define principal city as Philadelphia, and all other counties as suburbs
#GEOID for Philly county = 42101
long_MacroPM <- long_pm %>%  #long format
  mutate(
    macro_level = case_when(
      substr(GEOID, 1,5) %in% c("42101","34007", "10003") ~ "city",
      TRUE ~ "suburb"
 )
  )  

PM_data <- PM_data %>% #wide format (for divergence)
  mutate(
    macro_level = case_when(
      substr(GEOID, 1,5) %in% c("42101","34007", "10003") ~ "city",
      TRUE ~ "suburb"
    )
  )

#4.1 Dissimilarity on a macro level - city vs. suburbs？
D_macro <-segregation::dissimilarity(
  long_MacroPM,
  group = "tenure",
  unit = "macro_level",
  weight = "n"
)


#4.2 Divergence on Macro Level
#We just need to introduce grouping by places.
Div_macro <- PM_data %>%
  group_by(macro_level) %>%
  summarise(
    divergence_macro = divergence(
      percent_renters,
      percent_owners,
      population = totalE,
      comparison = comparison_tenure,
      summed = TRUE
    )
  )

#4.3 Theil's H Index on Macro Level
#Use rsegregation package; Tried the segrgegation package which still involves
#a lot manual calculations (very easy to get messy).
macro_summary <- PM_data %>% #Re-group by macro level
  group_by(macro_level) %>%
  summarise(
    renters = sum(renterE),
    owners = sum(ownerE),
    total = sum(totalE),
    .groups = "drop"
  ) %>% #percent columns deleted after the drop operation; manually adding them back
  mutate(
    percent_rent = renters / total,
    percent_own = owners / total
  )

H_between_macro <- rsegregation::entropy(
  macro_summary$percent_rent,
  macro_summary$percent_own,
  population = macro_summary$total,
  comparison = comparison_tenure,
  entropy_type = "information_theory",
  logBase = 2,
  summed = TRUE
)


#Overall Summary
suppressMessages({
  suppressWarnings({
    cat("\nMSA-level summary:\n")
    print(summary)
    
    cat("\nMacro-level Dissimilarity:\n")
    print(D_macro)
    
    cat("\nMacro-level Divergence:\n")
    print(Div_macro)
    
    cat("\nMacro-level Theil's H:\n")
    print(H_between_macro)
  }
  )
})

