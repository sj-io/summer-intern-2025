library(tidyverse)
library(tidycensus)
library(rsegregation) 

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
  filter(!is.na(pct_poverty)) %>%
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
      signed_di = (ai_over_AT - bi_over_BT)*100, #(multiplied by 100 to work with larger numbers)
      abs_di = abs(signed_di)
    ) %>%
    ungroup()
}
tract_di <- calculate_di(cleaned_set)
#2. convert into percentile cross MSA
add_di_percentile <- function(df) {
  df %>%
    group_by(year) %>%    mutate(
      di_percentile = percent_rank(abs_di) * 100
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

tract_di <- add_di_percentile(tract_di) %>% add_signed_percentile() %>% add_di_scaled()

#Modelling
#1. Prediction based
#1.1 Cleaning outcomes
#----
#####Juxtaposing the 4 decades - creating lagged outcome
# 1. Extract baseline percentiles/ signed_di
base80 <- tract_di %>% filter(year == 1980) %>% select(GEOID, di80 = signed_di)
base90  <- tract_di %>% filter(year == 1990) %>% select(GEOID, di90  = signed_di)
base00  <- tract_di %>% filter(year == 2000) %>% select(GEOID, di00  = signed_di)
base10  <- tract_di %>% filter(year == 2010) %>% select(GEOID, di10  = signed_di)
base19  <- tract_di %>% filter(year == 2019) %>% select(GEOID, di19  = signed_di)

# 2. Compute treatment flags for each decade window
t80_90 <- tract_di %>%
  filter(year == 1990) %>%   # LIHTC started in 1987
  group_by(GEOID) %>%
  summarise(t80_90 = as.integer(any(lihtc_placed == 1)), .groups="drop")

t91_00  <- tract_di %>% 
  filter(year == 2000) %>% 
  group_by(GEOID) %>% 
  summarise(t91_00 = as.integer(any(lihtc_placed == 1)), .groups="drop")

t01_10  <- tract_di %>% 
  filter(year == 2010) %>% 
  group_by(GEOID) %>% 
  summarise(t01_10 = as.integer(any(lihtc_placed == 1)), .groups="drop")

t11_19  <- tract_di %>% 
  filter(year == 2019) %>% 
  group_by(GEOID) %>% 
  summarise(t11_19 = as.integer(any(lihtc_placed == 1)), .groups="drop")

# 3. Stitch together one long data frame
plot_data_2.1 <- bind_rows(
  base80 %>% 
    left_join(t80_90, by = "GEOID") %>% 
    mutate(period = "1980", treated = t80_90, di = di80),
  
  base90 %>% 
    left_join(t91_00, by = "GEOID") %>% 
    mutate(period = "1990", treated = t91_00, di = di90),
  
  base00 %>% 
    left_join(t01_10, by = "GEOID") %>% 
    mutate(period = "2000", treated = t01_10, di = di00),
  
  base10 %>% 
    left_join(t11_19, by = "GEOID") %>% 
    mutate(period = "2010", treated = t11_19, di = di10),
) %>%
  mutate(
    treated = factor(treated, levels = c(0,1), labels = c("No LIHTC", "LIHTC"))
  ) %>%
  group_by(GEOID) %>%
  filter(n()==4) %>% 
  ungroup()




#2. Preparing predictors
#1.2 Preparing predictors
#----
predictors <- tract_di %>% 
  mutate(
    log_hinc = log(hinc),
    pct_poverty = pct_poverty*100,
  ) %>%
  mutate(race_group = case_when(
    pct_nhwht > 0.5 ~ "White",
    pct_nhblk > 0.5 ~ "Black",
    pct_hisp  > 0.5 ~ "Hispanic",
    pct_asian > 0.5 ~ "Asian",
    TRUE            ~ "Other"
  )) %>% mutate(
    year = as.character(year)
  ) %>% 
  select(GEOID, year, pct_poverty, log_hinc, race_group, cbsa13)

#Create 1980 racial group as a time-invariant variable
racial_majority_1980 <- predictors %>%
  filter(year == "1980") %>%
  select(GEOID, racial_majority_1980 = race_group)



#3. Merge two sets
merged_df <- plot_data_2.1%>%
  left_join(predictors, by = c("GEOID", "period" = "year"))
merged_df <- merged_df %>%
  mutate(decade = as.factor(period))
merged_df <- merged_df %>%
  mutate(lihtc_dummy = as.numeric(treated == "LIHTC"))

#1.3 Modeling
library(fixest)
merged_df$race_group <- relevel(factor(merged_df$race_group), ref = "White")
lpm_model <- feols(
  lihtc_dummy ~ di + pct_poverty + log_hinc + race_group + decade | cbsa13,
  data = merged_df
)
summary(lpm_model)

logit_model <- feglm(
  lihtc_dummy ~ di + pct_poverty + log_hinc + race_group + decade | cbsa13,
  data = merged_df,
  family = "binomial"
)
summary(logit_model)

dict_Q1 = c(
  "di" = "Signed DI",
  "pct_poverty" = "Poverty Rate (%)",
  "log_hinc" = "Log(Median Household Income)",
  "race_groupBlack" = "Majority: Black",
  "race_groupHispanic" = "Majority: Hispanic",
  "race_groupAsian" = "Majority: Asian",
  "race_groupOther" = "Majority: Others/No Majority",
  "decade1990" = "Decade: 1990",
  "decade2000" = "Decade: 2000",
  "decade2010" = "Decade: 2010"
)

etable(lpm_model, logit_model, dict= dict, tex = TRUE,file = "output/Q1RegressionModels.tex")

###Visualization
library(marginaleffects)
library(modelsummary)
coefs <- modelsummary::tidy(logit_model, conf.int = TRUE) %>% mutate(term_clean = dict_Q1[term], term_clean = factor(term_clean, levels = rev(dict_Q1))) 
ggplot(coefs, aes(x = estimate, y = term_clean)) +
  geom_point(size = 2.2) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.15) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  theme_bw(base_size = 12) +
  labs(
    x = "Log-Odds Coefficient Estimate (±95% CI)",
    y = NULL,
    title = "Logistic Model: Signed DI as a Predictor of LIHTC Placement", 
    caption = "Reference groups: Majority White; Decade = 1980"
  ) +  theme(
    axis.text.y = element_text(hjust = 0, face = "bold"),  #flatten look
    plot.title = element_text(face = "bold", size = 13, hjust = 1.35),
    plot.caption = element_text(size = 10, hjust = 0, color = "gray30"),
    panel.grid.major.y = element_blank()
  )
ggsave("output/Q1- 100MSA-PFERM.png", dpi = 300)


#Case Studies
NYC_df <- merged_df %>% 
  filter(cbsa13 == "35620") 
Philly_df <- merged_df %>% 
  filter(cbsa13 == "37980")
Houston_df <- merged_df %>% 
  filter(cbsa13 == "26420")
Chicago_df <- merged_df %>% 
  filter(cbsa13 == "16980")
NYC <- feglm(
  lihtc_dummy ~ di + pct_poverty + log_hinc + race_group + decade,
  data = NYC_df,
  family = "binomial"
)
Philadelphia <- feglm(
  lihtc_dummy ~ di + pct_poverty + log_hinc + race_group + decade,
  data = Philly_df,
  family = "binomial"
)
Houston <- feglm(
  lihtc_dummy ~ di + pct_poverty + log_hinc + race_group + decade,
  data = Houston_df,
  family = "binomial"
)
Chicago <- feglm(
  lihtc_dummy ~ di + pct_poverty + log_hinc + race_group + decade,
  data = Chicago_df,
  family = "binomial"
)
etable(NYC, Philadelphia, Chicago, Houston, dict = dict_Q1, tex = TRUE, file = "output/Q1PERM-CaseStudies.tex")

#VISUALIZATION

extract_di <- function(model, city_name) {
  modelsummary::tidy(model, conf.int = TRUE) %>%
    filter(term == "di") %>%
    mutate(city = city_name)
}

di_coefs <- bind_rows(
  extract_di(NYC, "New York"),
  extract_di(Philadelphia, "Philadelphia"),
  extract_di(Houston, "Houston"),
  extract_di(Chicago, "Chicago")
)


#2 Outcome focused
#----
library(fixest)
model_data <- tract_di %>% select(GEOID, year, signed_di, type, lihtc_placed) %>%
  mutate(year = as.numeric(year))
#2.1 constructing lithc dummy
model_data <- model_data %>%
  group_by(GEOID) %>%
  arrange(year) %>%
  mutate(lihtc_ever = as.numeric(cumsum(lihtc_placed == 1) > 0)) %>% #1.1 ever placed
  mutate(
    prev_lihtc = lag(cumsum(lihtc_placed == 1), default = 0),
    lihtc_first_time = ifelse(prev_lihtc == 0 & lihtc_placed == 1, 1, 0) #1.2 placed this decade AND previous decade
  ) %>%
  ungroup()


predictors <- predictors %>%
  mutate(year = as.numeric(year))
#join 2 datasets
model_data <- model_data %>% left_join(predictors, by = c("GEOID", "year"))

model_data <- model_data %>%
  left_join(racial_majority_1980, by = "GEOID")
#modelling
model_data$racial_majority_1980 <- relevel(factor(model_data$racial_majority_1980), ref = "White")
model_data$racial_group <- relevel(factor(model_data$racial_group), ref = "White")
outcome_model_1 <- feols(
  signed_di ~ lihtc_ever * racial_majority_1980 + pct_poverty + log_hinc  | GEOID,
  data = model_data
)
#################
####!!!!Use this alternative if you want overall patterns not separated by racial majority!!!!####
#outcome_model_1 <- feols(
#  signed_di ~ lihtc_ever + pct_poverty + log_hinc +racial_group | GEOID,
#  data = model_data
#)
############

outcome_model_2  <- feols(
  signed_di ~ lihtc_placed + pct_poverty + log_hinc + race_group | GEOID,
  data = model_data)

#Dont use this (bad logic)
#outcome_model_3 <- feols(
#  signed_di ~ lihtc_first_time + pct_poverty + log_hinc + race_group | GEOID,
#  data = model_data)

dict = c("lihtc_ever" = "LIHTC Placement (Cumulative)", 
         "signed_di" = "Signed Di",
         "race_groupBlack" = "Majority: Black",
         "race_groupHispanic" = "Majority: Hispanic",
         "race_groupAsian" = "Majority: Asian (Ref = White)",
         "race_groupOther" = "Majority: Others/No Majority",
         "log_hinc" = "Log(Median Household Income)",
         "pct_poverty" = "Poverty Rate (%)",
         "lihtc_placed" = "New LIHTC Placement During Decade",
         "lihtc_ever:racial_majority_1980Asian" = "LIHTC x Asian",
         "lihtc_ever:racial_majority_1980Black" = "LIHTC x Black",
         "lihtc_ever:racial_majority_1980Hispanic" = "LIHTC x Hispanic",
         "lihtc_ever:racial_majority_1980Other" = "LIHTC x Other")

etable(outcome_model_1, outcome_model_2, dict = dict, tex = TRUE, file ="output/Q2PFERM-AllMSAs.tex")
#Visualization 
library(modelsummary)
coefs_Q2 <- broom::tidy(outcome_model_1, conf.int = TRUE) %>%
  filter(term %in% names(dict)) %>%
  mutate(
    term_clean = dict[term],
    term_clean = factor(term_clean, levels = rev(dict))  # reverse for Signed DI on top
  )

ggplot(coefs_Q2, aes(x = estimate, y = term_clean)) +
  geom_point(size = 2.2) +
  geom_text(aes(label = round(estimate, 3), hjust = 0.5), 
            vjust = -0.7, size = 3.5) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.15) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  theme_bw(base_size = 12) +
  labs(
    x = "OLS Coefficient Estimates (±95% CI)",
    y = NULL,
    title = "Effect of LIHTC Placement on Signed DI", caption = "Reference groups: Majority White; Decade = 1980"
  ) +  theme(
    axis.text.y = element_text(hjust = 0, face = "bold"),     
    plot.title = element_text(face = "bold", size = 13),
    plot.caption = element_text(size = 10, hjust = 0, color = "gray30"),
    panel.grid.major.y = element_blank()
  )
ggsave("output/Q2-100MSA-PFERM-followup.png", dpi = 300)
#2.1Case Studies
#----
NYC_df <- model_data %>% 
  filter(cbsa13 == "35620") 
Philly_df <- model_data %>% 
  filter(cbsa13 == "37980")
Houston_df <- model_data %>% 
  filter(cbsa13 == "26420")
Chicago_df <- model_data %>% 
  filter(cbsa13 == "16980")
#1. By cumulative
NYC <- feols(
  signed_di ~ lihtc_ever + pct_poverty + log_hinc + race_group | GEOID,
  data = NYC_df
)
Philadelphia <- feols(
  signed_di ~ lihtc_ever + pct_poverty + log_hinc + race_group | GEOID,
  data = Philly_df
)
Houston <- feols(
  signed_di ~ lihtc_ever + pct_poverty + log_hinc + race_group | GEOID,
  data = Houston_df
)
Chicago <- feols(
  signed_di ~ lihtc_ever + pct_poverty + log_hinc + race_group | GEOID,
  data = Chicago_df
)
etable(NYC, Philadelphia, Houston, Chicago, dict = dict, tex = TRUE, file = "output/Q2PERM-CaseStudies.tex")

#2.2 Urban vs Suburban
urban_tracts <- model_data %>% filter(type == "Urban")
suburban_tracts <- model_data %>% filter(type == "Suburban")

#Option1
Urban_1 <- feols(signed_di ~ lihtc_ever + pct_poverty + log_hinc + race_group | GEOID, data = urban_tracts)
Suburban_1 <- feols(signed_di ~ lihtc_ever + pct_poverty + log_hinc + race_group | GEOID, data = suburban_tracts)
etable(Urban_1, Suburban_1, dict = dict, tex = TRUE, file = "output/Q2.2.1UrbanSuburban1.tex")
#Option 2
Urban_2 <- feols(signed_di ~ lihtc_placed + pct_poverty + log_hinc + race_group | GEOID, data = urban_tracts)
Suburban_2 <- feols(signed_di ~ lihtc_placed + pct_poverty + log_hinc + race_group | GEOID, data = suburban_tracts)
etable(Urban_2, Suburban_2, dict = dict, tex = TRUE, file = "output/Q2.2.2UrbanSuburban2.tex")

#2.3 By racial majority
racial_1980 <- model_data %>%
  filter(year == 1980) %>%
  select(GEOID, race_group_1980 = race_group)
model_data_23 <- model_data %>%
  left_join(racial_1980, by = "GEOID")
white_t <- model_data_23 %>% filter(race_group_1980 == "White")
black_t <- model_data_23 %>% filter(race_group_1980 == "Black")
hisp_t <- model_data_23 %>% filter(race_group_1980 == "Hispanic")
asian_t <- model_data_23 %>% filter(race_group_1980 == "Asian")
other_t <- model_data_23 %>% filter(race_group_1980 == "Other")

#run separate models
Majority_White <- feols(signed_di ~ lihtc_placed + pct_poverty + log_hinc | GEOID, data = white_t)
Majority_Black <- feols(signed_di ~ lihtc_placed + pct_poverty + log_hinc | GEOID, data = black_t)
Majority_Asian <- feols(signed_di ~ lihtc_placed + pct_poverty + log_hinc | GEOID, data = asian_t)
Majority_Hispanic <- feols(signed_di ~ lihtc_placed + pct_poverty + log_hinc | GEOID, data = hisp_t)
Majority_Other <- feols(signed_di ~ lihtc_placed + pct_poverty + log_hinc | GEOID, data = other_t)

etable(Majority_White, Majority_Black, Majority_Asian, Majority_Hispanic, Majority_Other, dict=dict,file = "output/Q2.3.2RaceGroups.tex")
