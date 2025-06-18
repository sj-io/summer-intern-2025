library(tidycensus)
library(tidyverse)
library(here)

#GOALS: On YEAR 2022
#1. total housing unit = 140.9 million
#2. 33.6% occupied by renters
#3. 29, 561 tracts are rental deserts (35.1% of all)
#4. 11, 313 tracts are extreme rental deserts (13.6% of all)

#1. Load variables to check what we might need
v_var <- load_variables(2022, "acs5")
#View(v_var)
#total housing unit = B25001_01
#total renter occupied = B25003_003
#total owner occupied = B25003_002
#vacant for rent = B25004_002
variables <-c(
  renters = "B25003_003",
  r_vacant = "B25004_002",
  total_housing = "B25001_001"
)

#2. Obtaining census data; need tract-level data for ALL states

all_states <-unique(fips_codes$state)[1:51]
data <- map_dfr(all_states, ~{
  get_acs(
    geography = "tract",
    variables = variables,
    year = 2022,
    survey = "acs5",
    state = .x,
    output = "wide",
  ) %>%
    mutate(state = .x
  ) %>% filter( ###Data Cleaning##
      !is.na(total_housingE),
      !is.na(rentersE),
      !is.na(r_vacantE),
      total_housingE > 0 #considering that some tracts report 0 for housing units, which cannot be used as a denominator for later calculation.
    ) %>% select(##only selecting the useful variables to reduce file size##
      GEOID, 
      state, 
      total_housingE,
      rentersE,
      r_vacantE
      ) 
})

write.csv(data, here("Airgood-Obrycki_data.rds"),row.names = FALSE) 
nation_data <- read.csv(here("Airgood-Obrycki_data.rds"))

#3. Computing National Characteristics
total_units <-sum(nation_data$total_housingE) #calculates total housing units
total_renterUnits = sum(nation_data$rentersE) + sum(nation_data$r_vacantE)
nation_data <- nation_data %>% mutate( #calculates rental share for every row
  rental_share = ((nation_data$r_vacantE+nation_data$rentersE)/nation_data$total_housingE)
) %>% 
  filter(!is.na(rental_share)
  )
RentalShare <- total_renterUnits/total_units #calculate nation-wide rental share

#4. Investigating rental deserts
total_tracts <-nrow(nation_data)
#Case 1: mild rental deserts
rental_deserts <- filter(nation_data, rental_share < 0.2)
percent_desert <-nrow(rental_deserts) %>% `/` (total_tracts)
total_RD <-nrow(rental_deserts)
#Case 2: extreme rental deserts
extreme_deserts <- filter(nation_data, rental_share < 0.1)
percent_extreme <-nrow(extreme_deserts) %>% `/` (total_tracts)
total_ED <-nrow(extreme_deserts)

#5. Data Presentation
summary_table <-tibble(
  categories = c(
    "Year",
    "Total Housing Units",
    "Percent of Total Rental Share",
    "Number of Rental Deserts",
    "Percent of Rental Deserts",
    "Number of Extreme Deserts",
    "Percent of Extreme Deserts"
  ),
  value = c(
    2022,
    total_units,
    RentalShare*100,
    total_RD,
    percent_desert*100,
    total_ED,
    percent_extreme*100
  )
)

print(summary_table)
