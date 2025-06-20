library(segregation)
library(tidycensus)
library(sf)
library(tidyverse)
Sys.getenv("CENSUS_API_KEY")  
 
philly_msa_fips <-"37980"

#Housing Tenure Variables
housing_variables <-c(
  total_occupied_housing= "B25003_001",
  owner_occupied = "B25003_002", 
  renter_occupied= "B25003_003" 
)
#Philly MSA area includes: PA-NJ-DE-MD.
#PA COUNTY FIPS
pa_counties <-c("017","029","045","091","101") #Pa: Bucks,Chester,Delaware, Montgomery, Philadelphia.  
#NJ COUNTY FIPS
nj_counties <-c("005","007","015","033") #NJ:Burlington, Camden,Gloucester, and Salem.
#DE COUNTY FIPS
de_counties <-c("003") #DE: New castle county 
#MD COUNTY FIPS
md_counties <-c("015") #MD: Cecil county.  

#Getting data for PA counties 
philly_data_pa <- get_acs(
  geography= "tract",
  variables = housing_variables,
  year= 2022, 
  survey = "acs5", 
  state= "PA", 
  county=pa_counties,
  output="wide", 
  geometry = FALSE
) 
  
#Getting data for NJ counties  
philly_data_nj <- get_acs(
  geography= "tract",
  variables = housing_variables,
  year= 2022, 
  survey = "acs5", 
  state= "NJ", 
  county=nj_counties,
  output="wide", 
  geometry = FALSE
) 

#Getting data for DE counties  
philly_data_de <- get_acs(
  geography= "tract",
  variables = housing_variables,
  year= 2022, 
  survey = "acs5", 
  state= "DE", 
  county=de_counties,
  output="wide", 
  geometry = FALSE
) 
#Getting data for MD counties  
philly_data_md <- get_acs(
  geography= "tract",
  variables = housing_variables,
  year= 2022, 
  survey = "acs5", 
  state= "MD", 
  county=md_counties,
  output="wide", 
  geometry = FALSE
) 

# Divergence Index (D)
#Here i am making the proportions for the whole msa 
#combining the different data frames I just made 

philly_msa_data_raw <-bind_rows (
  philly_data_pa,
  philly_data_nj,
   philly_data_de,
  philly_data_md
)

phl_tenure_data_clean <- philly_msa_data_raw %>% 
  select(GEOID,NAME,
    total_occupied= total_occupied_housingE,
    owner_occupied= owner_occupiedE,
    renter_occupied= renter_occupiedE 
    )%>%
  filter(total_occupied >0) %>%
 mutate( 
   owner_occupied=replace_na(owner_occupied, 0), 
   renter_occupied=replace_na(renter_occupied, 0) 
 ) 

#Making the data in the long format for the different indices
philly_msa_data <- phl_tenure_data_clean %>% 
  pivot_longer(
    cols= c(owner_occupied, renter_occupied),
    names_to= "tenure_group",
    values_to= "count"
  ) 

#MSA levels for divergence
msa_total_occupied <-sum(phl_tenure_data_clean$total_occupied,na.rm=TRUE)
msa_owner_total <-sum(phl_tenure_data_clean$owner_occupied,na.rm=TRUE)
msa_renter_total <-sum(phl_tenure_data_clean$renter_occupied,na.rm=TRUE)

msa_owner_prop <- ifelse(msa_total_occupied >0, msa_owner_total /msa_total_occupied,0) 
msa_renter_prop <-ifelse(msa_total_occupied >0, msa_renter_total /msa_total_occupied,0) 

#Here i am actually finding the value 
# p_owner_i = owner_occupied/total occupied 
# p_renter_i = renter_occupied/total occupied 

divergence_index <- phl_tenure_data_clean %>% 
  mutate( 
    term_owner= ifelse(owner_occupied > 0 & msa_owner_prop>0,
                       (owner_occupied/msa_owner_prop)* log(owner_occupied/(total_occupied *msa_owner_prop)), 0), 
    term_renter= ifelse(renter_occupied > 0 & msa_renter_prop >0, 
                        (renter_occupied/msa_renter_prop) * log(renter_occupied/(total_occupied *msa_renter_prop)),0),
    
    tract_contribution_divergence =  ifelse(msa_total_occupied >0, 
                                            (total_occupied/msa_total_occupied) *(term_owner + term_renter),0)
  ) %>% 
  summarise(Divergence_Index= sum(tract_contribution_divergence, na.rm=TRUE)) %>%
  pull(Divergence_Index)
cat("\n---Divergence Index (D) ---\n")
cat("Divergence Index (D):", divergence_index, "\n")
cat("measures the divergence in tenure between tracts and total MSA and higher value means more segregation.\n\n")

#Dissimilarity Index (D)
  dissimilarity_index <-philly_msa_data %>%
    filter(tenure_group %in% c("owner_occupied","renter_occupied")) %>% 
    dissimilarity ( 
    group= "tenure_group", 
    unit= "GEOID", 
    weight= "count"
  )
cat("\n---Dissimilarity Index (D) ---\n")
print(dissimilarity_index)
cat("The percentage value of one group of either renters or owners that would need to move to have an even spread in the MSA .\n\n")
cat ("\n")

#Isolation Index (P*)

isolation_indices_all <- philly_msa_data  %>% 
   isolation (
  group= "tenure_group", 
  unit= "GEOID", 
  weight= "count"
) 

print(isolation_indices_all) 
group_col_name <- NULL
for(col_name in names(isolation_indices_all)) {
  if (is.character(isolation_indices_all[[col_name]]) &&
      any(c("owner_occupied","renter_occupied")%in% isolation_indices_all[[col_name]])) {
    group_col_name <-col_name
    break
  }
}
if(is.null(group_col_name)) {
  stop("cannot find the column")
  } else { 
   cat(paste0("identified isolation group column as :'", group_col_name, "'\n"))
}

estimate_col_name_iso <- NULL
for(col_name in names(isolation_indices_all)) {
  if(is.numeric(isolation_indices_all[[col_name]]) &&
    !grepl("se|upper|lower|ci|pop", tolower(col_name)) ) {
    estimate_col_name_iso <-col_name
    break
  }
}
if(is.null(estimate_col_name_iso)) {
  stop("cannot find the estimate")
} else { 
  cat(paste0("identified isolation estimate column as :'", estimate_col_name_iso, "'\n"))
}

isolation_index_owner <- isolation_indices_all %>%
  filter(.data[[group_col_name]] == "owner_occupied")%>% 
  pull(.data[[estimate_col_name_iso]])

isolation_index_renter <- isolation_indices_all %>%
  filter(.data[[group_col_name]] == "renter_occupied")%>% 
  pull(.data[[estimate_col_name_iso]])

cat("\n---Isolation Index (P*)---\n")
cat("Isolation Index (P*) for Owners:", isolation_index_owner , "\n")
cat("Isolation Index (P*) for Renters:", isolation_index_renter, "\n")
cat ("\n")

cat("The probability that someone in a renter/owner census tract is the same.\n\n")


# Binary Information Index (H) or Theil's Entropy Index 

binary_info_index <- philly_msa_data %>% 
  mutual_total( 
    group = "tenure_group",
    unit= "GEOID",
    weight="count"
  ) 
cat("\n--- Binary Information Index (H) or Theil's Entropy Index ---\n")
cat ("\n")
print(binary_info_index) 
# measures overall residential segregation and values closer to 1 means more segregation")
  
