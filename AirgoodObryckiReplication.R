library(tidyverse)
library(dplyr)

acs_vars <- c(
  
  total_units="B25001_001", 
  renter_occupied= "B25003_003",
  vacant_for_rent= "B25004_003"
)

all_states <-c(state.abb, "DC", "PR")
housing_data <- NULL
for(st in all_states) {
message(paste("Downloading data for", st, ".."))

state_data <-tryCatch(
  {
    get_acs(
      geography = "tract",
      variables= acs_vars, 
      year=2022,
      survey="acs5", 
      output = "wide", 
      state= st
  )
  }, 
  error= function(e) {
    message(paste("Error Downloading data for", st, ":", e$message))
    return(NULL)
  })
if(!is.null(state_data)) {
  housing_data <-bind_rows(housing_data,state_data)

}
}
if(is.null(housing_data)|| nrow(housing_data)==0 ){
  stop("Housing data was not downloaded.")
}
housing_data <-housing_data %>%
  mutate( 
    rental_share=(renter_occupiedE + vacant_for_rentE)/total_unitsE,
rental_desert= rental_share < 0.2, 
extreme_rental_desert= rental_share < 0.1
) 
summary_replication <- housing_data %>%
  summarise(
    average_rental_share_across_tracts = mean(rental_share, na.rm=TRUE) *100, 
    total_housing_units_in_us=sum(total_unitsE, na.rm=TRUE),
    number_rental_deserts= sum(rental_desert, na.rm=TRUE), 
    total_nation_tracts = n(), 
    perecentage_rental_deserts= (number_rental_deserts/total_nation_tracts)*100, 
    number_extreme_rental_deserts= sum(extreme_rental_desert, na.rm=TRUE),
    percentage_extreme_rental_deserts= (number_extreme_rental_deserts/total_nation_tracts) *100
  )
# This lets see if my replication is like the findings in the article testing if On average, 33.6% of the housing in tracts is either occupied by a renter or vacant for rent, out of a total 140.9 million housing units in 2022. In 29,251 tracts, less than 20% of the housing stock is either renter- occupied or vacant for rent. These rental deserts account for 35.1% of the more than 84,000 tracts nationally but cover about two thirds of the country’s land area. In 11,313 tracts (13.6%), are what we consider extreme rental deserts in which less than 10% of the stock is occupied by renters or vacant for rent.

print(summary_replication, width=Inf)

# My results are slightly different as I got 31.7 percent instead of 33.6. 
# I got 142506742 and the article had 140.9 million housing units in 2022 
# I got 31,520 for the number of rental deserts, 85,396 total national tracts, and 36.9% of rental deserts 
# I got 12,472 for extreme rental deserts and then 14.6 for the percentage of extreme rental deserts. 

