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
#Ran a quick check for overlap between LTDB records with NA cbsa13 value and LIHTC. 
#Proved there are overlap, so we don't want to clean them now. 
subset_na <- ltdb %>%
  filter(is.na(cbsa13)) %>%
  rename(tractID = GEOID)

overlap <- inner_join(subset_na, lihtc, by = "tractID")



#Joining two datasets
fullset <- ltdb %>% 
  left_join(lihtc_summary,
            by = c("GEOID" = "tractID",
                   "year" = "year")) %>%
  mutate(
    lihtc_placed   = if_else(total_units == 0 | is.na(total_units), 0, 1)
  )



selected_set <- fullset %>%
  dplyr::select(
    GEOID,cbsa13,year, owner, renter, occhh, pct_rent, pct_own, lihtc_placed, total_units, low_income_units
  ) %>% filter(!is.na(cbsa13)) %>%
  filter(!(year %in% c(2012, 2020)))


##This is for checking any duplicated rows
sorted_set <- selected_set %>%
  group_by(GEOID, year) %>%
  filter(n() > 1) %>%
  arrange(GEOID, year)

### Labeling total units
selected_set <- selected_set %>% mutate(total_units = (ifelse(lihtc_placed == 0, 0, total_units)))
###Count how many rows have negative values
negative <- selected_set %>% filter(owner < 0 | renter < 0 | occhh < 0)
selected_set <-selected_set %>% anti_join(negative, by = "GEOID") 

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
tract_di <- calculate_di(selected_set)
#2. convert into percentile within each MSA
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


###Data visualization
#1. Logic/Predictor - Did LIHTC units go to already segregated areas?
di_1990 <- tract_di %>%
  filter(year == 1990) %>%
  dplyr::select(cbsa13, GEOID, di_percentile_1990 = di_percentile)

lihtc_2000 <- tract_di %>%
  filter(year == 2000) %>%
  dplyr::select(GEOID, LIHTC_1991_2000 = lihtc_placed)

plot_data <- di_1990 %>%
  left_join(lihtc_2000, by = "GEOID") %>%
  mutate(LIHTC_1991_2000 = ifelse(is.na(LIHTC_1991_2000), 0, LIHTC_1991_2000))

ggplot(plot_data, aes(x = factor(LIHTC_1991_2000, labels = c("No LIHTC", "LIHTC")), 
                      y = di_percentile_1990)) +
  geom_boxplot(fill = "plum") +
  labs(
    x = "LIHTC Placement (1991–2000)",
    y = "1990 Tenure Dissimilarity Percentile (within MSA)",
    title = "Did LIHTC Go to Already Segregated Areas?"
  ) + 
  theme_minimal()
###Signed version
di_1990 <- tract_di %>%
  filter(year == 1990) %>%
  select(GEOID, cbsa13, signed_di_1990 = signed_di) %>%
  distinct(GEOID, .keep_all = TRUE)

lihtc_2000 <- tract_di %>%
  filter(year == 2000) %>%
  select(GEOID, LIHTC_1991_2000 = lihtc_placed) %>%
  distinct(GEOID, .keep_all = TRUE)

plot_data <- di_1990 %>%
  left_join(lihtc_2000, by = "GEOID") %>%
  mutate(LIHTC_1991_2000 = ifelse(is.na(LIHTC_1991_2000), 0, LIHTC_1991_2000))

ggplot(plot_data, aes(x = factor(LIHTC_1991_2000, labels = c("No LIHTC", "LIHTC")), 
                      y = signed_di_1990)) +
  geom_violin(fill = "plum", alpha = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  labs(
    x = "LIHTC Placement (1991–2000)",
    y = "Signed Tenure Dissimilarity in 1990",
    title = "Were LIHTC Units Placed in Renter- or Owner-Heavy Tracts?"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "#252a32") +
  theme_bw() + coord_cartesian(ylim = c(-0.01, 0.01)) +
  theme(plot.title = element_text(face = "bold"))

#2. Treatment focus
lihtc_2000 <- tract_di %>%
  filter(year == 2000) %>%
  select(GEOID, LIHTC_1990_2000 = lihtc_placed) %>%
  distinct(GEOID, .keep_all = TRUE)

tract_wide <- tract_di %>%
  filter(year %in% c(1990, 2000)) %>%
  select(GEOID, cbsa13, year, signed_di, di_percentile) %>%
  pivot_wider(
    names_from = year,
    values_from = c(signed_di, di_percentile)
  ) %>%
  drop_na(signed_di_1990, signed_di_2000)

plot_data2 <- tract_wide %>%
  left_join(lihtc_2000, by = "GEOID") %>%
  mutate(
    LIHTC_1990_2000  = ifelse(is.na(LIHTC_1990_2000 ), 0, LIHTC_1990_2000),
    delta_signed_di = signed_di_2000 - signed_di_1990,
    delta_percentile = di_percentile_2000 - di_percentile_1990
  )

#2.1 by signed di
ggplot(plot_data2, aes(x = factor(LIHTC_1990_2000, labels = c("No LIHTC", "LIHTC")),
                      y = delta_signed_di)) +
  geom_boxplot(fill = "lightblue", width = 0.7, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Did Tenure Segregation Change After LIHTC Placement?",
    x = "LIHTC Placement (1991–2000)",
    y = "Change in Signed Tenure Dissimilarity (2000 - 1990)"
  ) +
  theme_bw()  +
  theme(plot.title = element_text(face = "bold")) + coord_cartesian(ylim = c(-0.0025, 0.0025))

#2.2 by percentile

ggplot(plot_data2, aes(x = factor(LIHTC_1990_2000, labels = c("No LIHTC", "LIHTC")), y = delta_percentile)) +
  geom_boxplot(fill = "skyblue", outlier.shape = NA, width = 0.4) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_cartesian(ylim = c(-10, 10)) +
  labs(
    title = "Did Relative Tenure Segregation Shift After LIHTC?",
    x = "LIHTC Placement (1990–2000)",
    y = "Change in Tenure Segregation Percentile (1990 - 2000)"
  ) +
  theme_minimal(base_size = 13)

#2.3 long term effects of early placements
seg_change <- tract_di %>%
  filter(year %in% c(1980, 2019)) %>%
  select(GEOID, year, di_percentile) %>%
  pivot_wider(names_from = year, values_from = di_percentile, names_prefix = "di_") %>%
  mutate(delta_percentile = di_2019 - di_1980)

treated_tracts <- tract_di %>%
  filter(year == 1990) %>%
  group_by(GEOID) %>%
  summarise(treated = as.integer(any(lihtc_placed == 1)))

analysis_data <- seg_change %>%
  left_join(treated_tracts, by = "GEOID") %>%
  mutate(treated = ifelse(is.na(treated), 0, treated))

ggplot(analysis_data, aes(x = factor(treated, labels = c("No LIHTC", "LIHTC 1981–1990")),
                          y = delta_percentile)) +
  geom_boxplot(fill = "skyblue", outlier.shape = NA, width = 0.4) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Long-Term Change in Tenure Segregation (1980–2019)",
    x = "Early LIHTC Placement (1981–1990)",
    y = "Change in Segregation Percentile"
  ) +
  theme_minimal()

t.test(delta_percentile ~ treated, data = analysis_data)


plot_data2 %>%
  group_by(treated = ifelse(LIHTC_2010_2019 > 0, "LIHTC", "No LIHTC")) %>%
  summarise(mean_change = mean(delta_percentile),
            se = sd(delta_percentile)/sqrt(n())) %>%
  ggplot(aes(x = treated, y = mean_change)) +
  geom_col(fill = "skyblue") +
  geom_errorbar(aes(ymin = mean_change - se, ymax = mean_change + se), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Mean Change in Segregation Percentile (2010–2019)",
    y = "Mean Δ Percentile",
    x = "LIHTC Placement"
  ) +
  theme_minimal()

#----
#####Juxtaposing the 4 decades
# 1. Extract baseline percentiles/ signed_di
base80 <- tract_di %>% filter(year == 1980) %>% select(GEOID, di80 = signed_percentile)
base90  <- tract_di %>% filter(year == 1990) %>% select(GEOID, di90  = signed_percentile)
base00  <- tract_di %>% filter(year == 2000) %>% select(GEOID, di00  = signed_percentile)
base10  <- tract_di %>% filter(year == 2010) %>% select(GEOID, di10  = signed_percentile)
base19  <- tract_di %>% filter(year == 2019) %>% select(GEOID, di19  = signed_percentile)

# 2. Compute treatment flags for each prior‐decade window
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
plot_data_all <- bind_rows(
  base80 %>% 
    left_join(t80_90, by="GEOID") %>% 
    mutate(period = "1980", treated = 0, di = di80),
  base90 %>% 
    left_join(t80_90, by="GEOID") %>% 
    mutate(period = "1990", treated = t80_90, di = di90),
  base00 %>% 
    left_join(t91_00, by="GEOID") %>% 
    mutate(period = "2000", treated = t91_00, di = di00),
  base10 %>% 
    left_join(t01_10, by="GEOID") %>% 
    mutate(period = "2010", treated = t01_10, di = di10),
  base19 %>% 
    left_join(t11_19, by="GEOID") %>% 
    mutate(period = "2019", treated = t11_19, di = di19)
) %>%
  mutate(
    treated = factor(treated, levels = c(0,1), labels = c("No LIHTC","LIHTC"))
  )


plot_data_all <- plot_data_all %>%
  group_by(GEOID) %>%
  filter(n()==5) %>%
  ungroup()

# 4. Plot them side by side, grouped by treatment
ggplot(plot_data_all, aes(x = period, y = di, fill = treated)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.65, outlier.shape = 4, color = "black") +
  scale_fill_brewer("LIHTC in prior decade", palette="Set3") +
  labs(
    title = "Cross-MSA Tenure Segregation Percentile by LIHTC Exposure",
    subtitle = "Comparing 1990, 2000, 2010, 2019\n(each panel shows whether tract got LIHTC in the prior decade)",
    x = "Census Year",
    y = "Signed Percentile",
    caption = "(+100 = Renter-heavy tracts, –100 = Owner-heavy tracts)"
  ) +
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = round(..y.., 1)),
    position = position_dodge(width = 0.75),
    vjust = -0.5, size = 3.3, color = "black"
  ) +
  theme_minimal(base_size=15) +
  theme(legend.position = "top", legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 16, margin = margin(b = 6)),
        plot.subtitle = element_text(size = 13, margin = margin(b = 10)),
        plot.caption = element_text(size=12, hjust = 0.5),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))
  )

ggsave("LIHTC_boxplot.png", width = 9, height = 6, dpi = 300)

#----
#4.2 Min_max scale plot
ggplot(tract_di, aes(x = factor(year), y = di_scaled, fill = factor(lihtc_placed))) +
  geom_boxplot(outlier.shape = 1, position = position_dodge(width = 0.9), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_brewer("LIHTC in prior decade", palette="Set3") +
  scale_y_continuous(limits = c(-100, 100)) +
  labs(
    title = "Min–Max Scaled Signed Tenure Segregation",
    subtitle = "Comparing 1990, 2000, 2010, 2019\n(each panel shows whether tract got LIHTC in the prior decade)",
    x = "Census Year",
    y = "Scaled DI (-100 = Owner-heavy, +100 = Renter-heavy)",
    fill = "LIHTC in Prior Decade"
  ) +
  theme_minimal(base_size = 14) +
  stat_summary(
  fun = median,
  geom = "text",
  aes(label = round(..y.., 1)),
  position = position_dodge(width = 0.75),
  vjust = -0.5, size = 3.5, color = "black" ) +
  theme(legend.position = "top", legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 16, margin = margin(b = 6)),
        plot.subtitle = element_text(size = 13, margin = margin(b = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10))
  )
