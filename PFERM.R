library(tidycensus)
library(rsegregation)
library(tidyverse)
#Import the two datasets
ltdb <- read.csv("updated_ltdb.csv") %>% mutate(
  GEOID = str_pad(as.character(GEOID),width = 11, pad = "0")
)
lihtc <- read.csv("Revised_100MSA_lihtc_cleaned.csv") %>% mutate(
  tractID = str_pad(as.character(tractID), width = 11, pad = "0")
) 

#Data cleaning
#1. relabel decadal bucket into 2019
lihtc <- lihtc %>% mutate(
   year = if_else(
   year == 2020,
    2019,
    year
    )
)

#2. Aggregate LIHTC data so that each year+GEOID pair has only one row
lihtc_summary <- lihtc %>%
  group_by(tractID, year) %>%
  summarise(
    total_units = sum(total_units, na.rm = TRUE),
    low_income_units = sum(low_income_units, na.rm = TRUE),
    lihtc_placed = as.integer(total_units > 0)
  )

#3. supplementing the missing 2010 hinc, %foreign & %poverty values with 2012 data
hinc_pov_2012 <- ltdb %>%
  filter(year == 2012) %>%
  dplyr::select(GEOID, hinc,pct_foreign, pct_poverty) %>%
  rename(hinc_2012 = hinc, pct_foreign_2012 = pct_foreign, pct_poverty_2012 = pct_poverty)
ltdb_2010_joined <- ltdb %>%
  filter(year == 2010) %>%
  left_join(hinc_pov_2012, by = "GEOID") %>% mutate(hinc = hinc_2012, pct_foreign = pct_foreign_2012, pct_poverty = pct_poverty_2012) %>%
  dplyr::select(-hinc_2012, -pct_poverty_2012, -pct_foreign_2012)
ltdb_other <- ltdb %>%
  filter(year != 2010)
ltdb <- bind_rows(ltdb_2010_joined, ltdb_other)
#Joining two datasets
fullset <- ltdb %>% 
  left_join(lihtc_summary,
            by = c("GEOID" = "tractID",
                   "year" = "year")) %>%
  mutate(
    lihtc_placed   = if_else(total_units == 0 | is.na(total_units), 0, 1)
  )



selected_set <- fullset %>%
 filter(!is.na(cbsa13)) %>% 
  filter(!is.na(hinc)) %>%
  filter(!(year %in% c(2020,2012)))
### Labeling total units
selected_set <- selected_set %>% mutate(total_units = (ifelse(lihtc_placed == 0, 0, total_units)))
###Count how many rows have negative valus
negative <- selected_set %>% filter(owner < 0 | renter < 0 | occhh < 0| hinc < 0 | npov < 0 | nhwht < 0 | nhblk < 0| hisp <0 | asian < 0 |ntv < 0 |fhh <0 | fb <0)
selected_set <-selected_set %>% anti_join(negative, by = "GEOID") 

##This is for checking any duplicated rows
sorted_set <- selected_set %>%
  group_by(GEOID, year) %>%
  filter(n() > 1) %>%
  arrange(GEOID, year)
#More data cleaning - checking errorneous data entry
problematic <- selected_set %>%
  filter(occhh > 0 & owner == 0 & renter == 0)

nrow(problematic)
bad_tracts <- unique(problematic$GEOID)
selected_set <- selected_set %>%
  filter(!(GEOID %in% bad_tracts))


#Creating fixed panel
cleaned_set <- selected_set %>%
  group_by(GEOID) %>%
  filter(n()== 5) %>%
  ungroup()
#Propose functions for calculating the dissimilarity index
#1. tract level composition d
calculate_di <- function(df) {
  df %>%
    group_by(cbsa13, year) %>%
    mutate(
      AT = sum(renter, na.rm = TRUE),
      BT = sum(owner, na.rm = TRUE),
      ai_over_AT = renter / AT,
      bi_over_BT = owner / BT,
      signed_di = ai_over_AT - bi_over_BT,
      abs_di = abs(signed_di)
    ) %>%
    ungroup()
}
tract_di <- calculate_di(cleaned_set)

#2. convert into percentile within each MSA
add_di_percentile <- function(df, di_col = "abs_di") {
  df %>%
    group_by(year) %>%
    mutate(
      di_percentile = percent_rank(.data[[di_col]]) * 100
    ) %>%
    ungroup()
}

add_di_scaled <- function(df) { #Cross msa comparison
  df %>%
    group_by(year) %>%
    mutate(
      max_abs_di = max(abs(signed_di), na.rm = TRUE),
      di_scaled = if_else(
        max_abs_di == 0,
        0,
        (signed_di / max_abs_di) * 100
      )
    ) %>%
    ungroup()
}
add_signed_percentile <- function(df) {
  df_neg <- df %>%
    filter(signed_di < 0) %>%
    group_by(year) %>%
    mutate(
      signed_percentile = -((rank(-signed_di, ties.method = "min") - 1) /
                              (n() - 1)) * 100
    ) %>%
    ungroup()
  
  df_pos <- df %>%
    filter(signed_di > 0) %>%
    group_by(year) %>%
    mutate(
      signed_percentile = ((rank(signed_di, ties.method = "min") - 1) /
                             (n() - 1)) * 100
    ) %>%
    ungroup()
  
  df_zero <- df %>%
    filter(signed_di == 0) %>%
    mutate(signed_percentile = 0)
  
  # Combine and return
  bind_rows(df_neg, df_zero, df_pos) %>%
    arrange(cbsa13, year, GEOID)
}


tract_di <- add_di_percentile(tract_di, di_col = "abs_di") %>% add_signed_percentile() %>% add_di_scaled()

#Case Studies
NYC <- tract_di %>% 
  filter(cbsa13 == "35620") 
Philly <- tract_di %>% 
  filter(cbsa13 == "37980")
Houston <- tract_di %>% 
  filter(cbsa13 == "26420")
Chicago <- tract_di %>% 
  filter(cbsa13 == "16980")

#Inspect the outcome variable with data visualization
library(patchwork)
library(sf)
library(fixest)
epsilon <- 1e-5 
mhv_histogram <- ggplot(tract_di, aes(x = signed_di)) + 
  geom_histogram(alpha = 0.5, fill = "navy", color = "navy",
                 bins = 100) + 
  theme_minimal() + 
  scale_x_continuous(labels = function(x) paste0(x, "%")) + 
  labs(x = "tract level contribution to dissimilarity value")

mhv_histogram

#Regression model
tract_di_check <- tract_di %>%
  mutate(
    y_log = log(abs_di + epsilon),
    sqrt_pov = sqrt(pct_poverty),
    white = pct_nhwht,
    sqrt_black = sqrt(pct_nhblk),
    sqrt_hisp = sqrt(pct_hisp),
    sqrt_asian = sqrt(pct_asian),
    sqrt_foreign = sqrt(pct_foreign),
    log_hinc = log(hinc)
  )  

model_data <- tract_di_check %>%
  filter(complete.cases(y_log, sqrt_pov, sqrt_asian, white, sqrt_hisp,
                         sqrt_black, sqrt_foreign, log_hinc)) %>%
  filter(is.finite(log_hinc)) 

model1<- feols(y_log ~ lihtc_placed + I(lihtc_placed * total_units) + sqrt_pov + sqrt_black +
                 sqrt_hisp + white +sqrt_asian + sqrt_foreign + log_hinc 
               | GEOID + year, cluster = ~GEOID, data = model_data)

model_data_descriptivestats <- model_data %>% dplyr::select(GEOID, year, cbsa13, pct_poverty,
                                                            pct_nhwht, pct_nhblk, pct_hisp, pct_asian, pct_foreign, 
                                                            total_units, lihtc_placed, abs_di, signed_di, di_percentile,
                                                            y_log, sqrt_pov, white, sqrt_black, sqrt_hisp, sqrt_asian,
                                                            sqrt_foreign, log_hinc)




NYCMSA <- model1
PhiladelphiaMSA <- model1
HoustonMSA <- model1
ChicagoMSA <- model1
model1
#Exporting data
library(webshot2)
library(flextable)
etable_out <- etable(model1,
               dict = c(
                 "y_log" = "Tract DI Component",
                 "sqrt_pov" = "Poverty Rate",
                 "sqrt_white" = "% Non-Hispanic White",
                 "sqrt_black" = "% Non-Hispanic Black",
                 "sqrt_hisp" = "% Hispanic",
                 "sqrt_asian" = "% Asian",
                 "white" = "% White",
                 "sqrt_foreign" = "% Foreign Born",
                 "log_hinc" = "Median Household Income",
                 "I(lihtc_placed * total_units)" = "LIHTC Units (only where placed)",
                 "lihtc_placed" = "LIHTC Placement (Yes/No)", tex = TRUE, file = "PFERM_tractdi.tex")
               )
  
etable_out

#######Introducing Theils's H##########
library(rsegregation)
# Filter MSAs with both urban and suburban tracts
msas_with_both <- tract_di %>%
  distinct(cbsa13, type) %>%
  group_by(cbsa13) %>%
  summarise(n_types = n_distinct(type), .groups = "drop") %>%
  filter(n_types == 2) %>%
  pull(cbsa13)

# Then use this list to filter the full dataset
tract_di_filtered <- tract_di_filtered %>%
  filter(cbsa13 %in% msas_with_both)

tract_summary <- tract_di_filtered %>%
  group_by(cbsa13, year, type) %>%
  summarise(
    renters = sum(renter, na.rm = TRUE),
    owners = sum(owner, na.rm = TRUE),
    total = sum(occhh, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_rent = renters / total,
    pct_own = owners / total
  ) 
msa_comparison <- tract_di_filtered %>%
  group_by(cbsa13, year) %>%
  summarise(
    renters = sum(renter, na.rm = TRUE),
    owners = sum(owner, na.rm = TRUE),
    total = sum(occhh, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    comp_rent = renters / total,
    comp_own = owners / total
  )

tract_summary <- tract_summary %>%
  left_join(msa_comparison, by = c("cbsa13", "year"))

H_by_msa_year <- tract_summary %>%
  group_split(cbsa13, year) %>%
  map_dfr(~{
    group_data <- .
    data.frame(
      cbsa = unique(group_data$cbsa13),
      year = unique(group_data$year),
      theil_H = rsegregation::entropy(
        group_data$pct_rent,
        group_data$pct_own,
        population = group_data$total.x,
        comparison = c(unique(group_data$comp_rent), unique(group_data$comp_own)),
        entropy_type = "information_theory",
        logBase = 2,
        summed = TRUE
      )
    )  })

#join back to lihtc and other covariates
macro_predictors <- tract_di %>%
  group_by(cbsa13, year) %>%
  summarise(
    total_units = sum(total_units, na.rm = TRUE),
    pct_poverty = weighted.mean(pct_poverty, occhh, na.rm = TRUE),
    pct_foreign = weighted.mean(pct_foreign, occhh, na.rm = TRUE),
    pct_nhblk = weighted.mean(pct_nhblk, occhh, na.rm = TRUE),
    pct_nhwht = weighted.mean(pct_nhwht, occhh, na.rm =TRUE),
    pct_hisp = weighted.mean(pct_hisp, occhh, na.rm =TRUE),
    pct_asian = weighted.mean(pct_asian, occhh, na.rm = TRUE),
    hinc = weighted.mean(hinc, occhh, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(cbsa = cbsa13)  #aggregate to msa level
panel_data <- H_by_msa_year %>%
  left_join(macro_predictors, by = c("cbsa", "year")) 


model_H <- feols(
  theil_H ~ total_units + pct_poverty + pct_foreign + hinc + pct_nhblk + pct_hisp + pct_asian + pct_nhwht | cbsa + year,
  data = panel_data,
  cluster = ~cbsa
)
summary(model_H)

#Formatting
etable_H <- etable(model_H,
                     dict = c(
                       "total_units" = "Total Number of LIHTC Units",
                       "pct_poverty" = "Poverty Rate",
                       "pct_nhwht" = "% Non-Hispanic White",
                       "pct_nhblk" = "% Non-Hispanic Black",
                       "pct_hisp" = "% Hispanic",
                       "pct_asian" = "% Asian",
                       "pct_foreign" = "% Foreign Born",
                       "hinc" = "Median Household Income",
                       "theil_H" = "Theil's H Value", tex = TRUE, file = "output/PFERM_H.tex"
))

etable_H



##useless
mhv_histogram <- ggplot(panel_data, aes(x = log(theil_H))) + 
  geom_histogram(alpha = 0.5, fill = "navy", color = "navy",
                 bins = 100) + 
  theme_minimal() + 
  scale_x_continuous(labels = function(x) paste0(x, "%")) + 
  labs(x = "tract level contribution to dissimilarity value")

mhv_histogram

###Data Visualization
library(patchwork)
library(carr)
plot_H <- ggplot(panel_data, aes(x = year, y = theil_H, color = cbsa, group = cbsa)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(
    title = "Theil's H Trend by MSA",
    x = "Year",
    y = "Theil's H (Tenure Segregation)",
    color = "MSA"
  ) +
  theme_minimal(base_size = 14)

plot_lihtc <- ggplot(panel_data, aes(x = year, y = total_units, fill = cbsa, group = cbsa)) +
  geom_col(position = "dodge") +
  labs(
    title = "LIHTC Units Placed by Year and MSA",
    x = "Year",
    y = "Total LIHTC Units",
    fill = "MSA"
  ) +
  theme_minimal(base_size = 14)

