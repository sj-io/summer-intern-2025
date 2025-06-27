#Goal: Clean LIHTC data
library(tidycensus)
library(tidyverse)
library(stringr)
LIHTC_raw <- read.csv("LIHTCPUB.csv")
LIHTC_raw$st2020 <- as.numeric(LIHTC_raw$st2020)
LIHTC_raw$cnty2020 <- as.numeric(LIHTC_raw$cnty2020)

fips_lookup <- read.csv("national_county.txt", sep = ",", header = FALSE, col.names = c("state_abbr", "state_fips", "county_fips", "county_name", "class_code"))

# Create full county FIPS code
fips_lookup$full_fips <- sprintf("%02d%03d", fips_lookup$state_fips, fips_lookup$county_fips)
LIHTC_raw$county_fips_5 <- ifelse(is.na(LIHTC_raw$st2020) | is.na(LIHTC_raw$cnty2020),
  NA, sprintf("%02d%03d", LIHTC_raw$st2020, LIHTC_raw$cnty2020))
  
glimpse(LIHTC_raw)
target_counties <- data.frame(
  state_abbr = c(rep("IL", 9), rep("IN", 4), "WI",rep("TX", 9), rep("NY", 12), rep("NJ", 12), "PA",rep("PA", 5), rep("NJ", 4), "DE", "MD"),
  county_name = c("Cook County", "DuPage County", "Grundy County", "Kendall County", "McHenry County", "Will County",
    "DeKalb County", "Kane County", "Lake County",
    "Jasper County", "Lake County", "Newton County", "Porter County",
    "Kenosha County",
    "Austin County", "Brazoria County", "Chambers County", "Fort Bend County",
    "Galveston County", "Harris County", "Liberty County", "Montgomery County", "Waller County",
    "Bronx County", "Dutchess County", "Kings County", "Nassau County", "New York County", "Orange County",
    "Putnam County", "Queens County", "Richmond County", "Rockland County", "Suffolk County", "Westchester County",
    "Bergen County", "Essex County", "Hudson County", "Hunterdon County", "Middlesex County", "Monmouth County",
    "Morris County", "Ocean County", "Passaic County", "Somerset County", "Sussex County", "Union County",
    "Pike County",
    "Bucks County", "Chester County", "Delaware County", "Montgomery County", "Philadelphia County",
    "Burlington County", "Camden County", "Gloucester County", "Salem County",
    "New Castle County", "Cecil County")
)
#A very tiresome way to filter out the counties that we want :).
target_fips <- merge(target_counties, fips_lookup, by = c("state_abbr", "county_name"))
lihtc_filtered <- LIHTC_raw %>% filter(county_fips_5%in% target_fips$full_fips)
view(lihtc_filtered)

#Now it's time to filter the variables we need
lihtc_cleaned <- lihtc_filtered %>% 
  filter(yr_pis <= 2020) %>%
  select(hud_id, project, proj_add,
         county_fips_5, fips2010, fips2020, yr_pis, n_units, li_units)


#Cleaning units with invalid 2010 tract ID via crosswalk
#I checked and all missing units have valid 2020 tractID
fips_crosswalk <- read.csv("crosswalk_2020_2010.csv")
fips_crosswalk$TRACT2020 <- substr(as.character(fips_crosswalk$BLKID2020), 1, 11)
fips_crosswalk <- fips_crosswalk %>%
  distinct(TRACT2020, .keep_all = TRUE)

missing_2010 <- lihtc_cleaned %>%
  filter(str_detect(fips2010, "X"))


missing_2010_joined <- missing_2010 %>%
  left_join(fips_crosswalk,
            by = c("fips2020" = "TRACT2020")) %>%
  mutate(fips2010 = TRTID2010)

#Overwrite the missing fips2010 slots 
lihtc_cleaned <- lihtc_cleaned %>%
  left_join(missing_2010_joined %>% select(hud_id, fips2010_new = fips2010), by = "hud_id") %>%
  mutate(
    fips2010 = ifelse(str_detect(fips2010, "X"), fips2010_new, fips2010)
  ) %>%
  select(-fips2010_new)

#Clean data with invalid li_units or unreasonable total_units (0 total when more than 0 li_units):
lihtc_cleaned$li_units <- ifelse(is.na(lihtc_cleaned$li_units), 0, lihtc_cleaned$li_units)
lihtc_cleaned <- lihtc_cleaned %>%
  mutate(n_units = ifelse(n_units == 0 & li_units > 0, li_units, n_units)) %>%
  filter(!is.na(n_units))

#Categorize LIHTC data into decennial buckets
lihtc_binned <- lihtc_cleaned %>%
  mutate(decade_bucket = case_when(
    yr_pis >= 1987 & yr_pis <= 1990 ~ "1990",
    yr_pis >= 1991 & yr_pis <= 2000 ~ "2000",
    yr_pis >= 2001 & yr_pis <= 2010 ~ "2010",
    yr_pis >= 2011 & yr_pis <= 2020 ~ "2020",
  ))


final_lihtc <-lihtc_binned %>% select(hud_id, fips2010,decade_bucket, n_units, li_units) %>%
  rename(tractID = fips2010,
         year = decade_bucket,
         total_units = n_units,
         low_income_units = li_units)


write.csv(final_lihtc,"lihtc_cleaned.csv",row.names = FALSE)
