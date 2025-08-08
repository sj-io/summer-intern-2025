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
#2. convert into percentile WIHTIN MSA
add_di_percentile <- function(df) {
  df %>%
    group_by(year,cbsa13) %>%    mutate(
      di_percentile = percent_rank(abs_di) * 100
    ) %>%
    ungroup()
}


add_signed_percentile <- function(df) {
  df_neg <- df %>%
    filter(signed_di < 0) %>%
    group_by(year,cbsa13) %>%
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



tract_di <- add_di_percentile(tract_di) %>% add_signed_percentile() 


###Data visualization
#----
#####Juxtaposing the 4 decades
# 1. Extract baseline percentiles/ signed_di
base80 <- tract_di %>% filter(year == 1980) %>% select(GEOID, di80 = signed_percentile,cbsa13)
base90  <- tract_di %>% filter(year == 1990) %>% select(GEOID, di90  = signed_percentile,cbsa13)
base00  <- tract_di %>% filter(year == 2000) %>% select(GEOID, di00  = signed_percentile,cbsa13)
base10  <- tract_di %>% filter(year == 2010) %>% select(GEOID, di10  = signed_percentile,cbsa13)
base19  <- tract_di %>% filter(year == 2019) %>% select(GEOID, di19  = signed_percentile,cbsa13)

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

#CASE STUDIES
plot_data_2.1 <- plot_data_2.1 %>%
  filter(cbsa13 %in% c("35620", "37980", "26420", "16980")) %>%
  mutate(
    case_study_label = recode(cbsa13,
                              "35620" = "New York-Newark-Jersey City",
                              "37980" = "Philadelphia-Camden-Wilmington",
                              "26420" = "Houston-The Woodlands-Sugar Land",
                              "16980" = "Chicago-Naperville-Elgin"
    )
  )

#Define Color Palette
city_colors <- list(
  "New York-Newark-Jersey City" = c("No LIHTC" = "#a6cee3", "LIHTC" = "#1f78b4"),
  "Philadelphia-Camden-Wilmington" = c("No LIHTC" = "#b2df8a", "LIHTC" = "#33a02c"),
  "Houston-The Woodlands-Sugar Land" = c("No LIHTC" = "#fb9a99", "LIHTC" = "#e31a1c"),
  "Chicago-Naperville-Elgin" = c("No LIHTC" = "#fdbf6f", "LIHTC" = "#ff7f00")
)

# 4. Plot them side by side, grouped by treatment
plot_city <- function(city_name) {
  df <- plot_data_2.1 %>% filter(case_study_label == city_name)
  ggplot(df, aes(x = period, y = di, fill = treated)) +
    geom_boxplot(position = position_dodge(width = 0.75), width = 0.65, outlier.shape = NA, color = "black") +
    scale_fill_manual("LIHTC in the Next Decade", values = city_colors[[city_name]]) +
    labs(
      title = paste0("Tenure Segregation in ", city_name),
      subtitle = "Comparing 1990, 2000, 2010, 2019\n(each panel shows whether tract got LIHTC in the following decade)",
      x = "Census Year",
      y = "Signed Percentile",
      fill = "LIHTC",
      caption = "(+100 = Renter-heavy tracts, –100 = Owner-heavy tracts)"
    ) +
    stat_summary(
      fun = median,
      geom = "text",
      aes(label = round(..y.., 1)),
      position = position_dodge(width = 0.75),
      vjust = -0.5, size = 3.3, color = "black"
    ) +     stat_summary(
      fun = median,
      geom = "text",
      aes(label = round(..y.., 1), group = treated),
      position = position_dodge(width = 0.75),
      vjust = -0.6,
      size = 3.3,
      color = "black"
    ) +     stat_summary(
      fun = median,
      geom = "line",
      aes(group = treated),
      position = position_dodge(width = 0.75),
      linewidth = 0.9,
      alpha = 0.4
    ) +
    theme_bw(base_size=14) +
    theme(
      legend.position = "top",legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12),
      plot.caption = element_text(size=12, hjust = 0.5),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.title.x = element_text(margin = margin(t = 10))
    )
}
ny_plot <- plot_city(  "New York-Newark-Jersey City")
#ggsave("output/Q1.1-NYC.png", width = 9, height = 6, dpi = 300)
#philly_plot <- plot_city("Philadelphia-Camden-Wilmington")
#ggsave("output/Q1.1-Philly.png", width = 9, height = 6, dpi = 300)
#houston_plot <- plot_city("Houston-The Woodlands-Sugar Land")
#ggsave("output/Q1.1-Houston.png", width = 9, height = 6, dpi = 300)
#chicago_plot <- plot_city("Chicago-Naperville-Elgin")
#ggsave("output/Q1.1-Chicago.png", width = 9, height = 6, dpi = 300)

#----
#### 4.2 Treatment focused
#defining delta changes
tract_wide <- tract_di %>%
      select(GEOID, cbsa13, year, signed_percentile,lihtc_placed) %>%
       pivot_wider(
            names_from = (year),
            values_from = c(signed_percentile, lihtc_placed),
            names_glue = "{.value}_{year}"
       )
       

tract_wide <-tract_wide %>% mutate(
  delta_8090 = signed_percentile_1990 - signed_percentile_1980,
  delta_9000 = signed_percentile_2000 -signed_percentile_1990,
  delta_0010 = signed_percentile_2010 - signed_percentile_2000,
  delta_1019 = signed_percentile_2019 - signed_percentile_2010
)

delta_long <- tract_wide %>%
  select(GEOID, cbsa13,
         delta_8090, lihtc_placed_1990,
         delta_9000, lihtc_placed_2000,
         delta_0010, lihtc_placed_2010,
         delta_1019, lihtc_placed_2019) %>%
  pivot_longer(
    cols = starts_with("delta_"),
    names_to = "period",
    names_prefix = "delta_",
    values_to = "delta_di"
  ) %>%
  mutate(
    treat = case_when(
      period == "8090" ~ lihtc_placed_1990,
      period == "9000" ~ lihtc_placed_2000,
      period == "0010" ~ lihtc_placed_2010,
      period == "1019" ~ lihtc_placed_2019
    ),
    treated = factor(treat, levels = c(0,1), labels = c("No LIHTC", "LIHTC")),
    period = case_when(
      period == "8090" ~ "1980–1990",
      period == "9000" ~ "1990–2000",
      period == "0010" ~ "2000–2010",
      period == "1019" ~ "2010–2019"
    )
  ) %>%
  filter(!is.na(delta_di))


###CASE STUDIES
delta_long_case <- delta_long %>%
  filter(cbsa13 %in% c("35620", "37980", "26420", "16980")) %>%
  mutate(
    case_study_label = recode(cbsa13,
                              "35620" = "New York-Newark-Jersey City",
                              "37980" = "Philadelphia-Camden-Wilmington",
                              "26420" = "Houston-The Woodlands-Sugar Land",
                              "16980" = "Chicago-Naperville-Elgin"
    ))

plot_delta_city <- function(city_name) {
  df <- delta_long_case %>% filter(case_study_label == city_name)
  ggplot(df, aes(x = period, y = delta_di, fill = treated)) +
    geom_boxplot(position = position_dodge(width = 0.7), width = 0.55, outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
    stat_summary(
      fun = median,
      geom = "text",
      aes(label = round(..y.., 1), group = treated),
      position = position_dodge(width = 0.75),
      vjust = -0.5,
      size = 3.3,
      color = "black"
    ) +
    
    scale_fill_manual(values = city_colors[[city_name]]) +
    coord_cartesian(ylim = c(-30, 30)) +
    scale_y_continuous(breaks = seq(-30, 30, by = 10)) +
    labs(
      title = paste0("Δ Signed DI in ", city_name),
      subtitle = "Change over each decade by LIHTC treatment",
      x = "Decade",
      y = "Change in DI (ΔDI)",
      fill = "LIHTC",
      linetype = "Trend"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top",
      legend.title = element_text(face = "bold")
    )
}

ny_delta     <- plot_delta_city("New York-Newark-Jersey City")
philly_delta <- plot_delta_city("Philadelphia-Camden-Wilmington")
houston_delta<- plot_delta_city("Houston-The Woodlands-Sugar Land")
chicago_delta<- plot_delta_city("Chicago-Naperville-Elgin")

#Exporting plots
ny_delta
ggsave("output/Q2.1-NYC.png", width = 9, height = 6, dpi = 300)
philly_delta
ggsave("output/Q2.1-Philly.png", width = 9, height = 6, dpi = 300)
houston_delta
ggsave("output/Q2.1-Houston.png", width = 9, height = 6, dpi = 300)
chicago_delta
ggsave("output/Q2.1-Chicago.png", width = 9, height = 6, dpi = 300)

#----
#4.1.2 Bar graphs
delta_summary <- delta_long %>%
  group_by(period, treated) %>%
  summarise(median_delta = median(delta_di, na.rm = TRUE), .groups = "drop")

ggplot(delta_summary, aes(x = period, y = median_delta, fill = treated)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values =c("No LIHTC" = "#bcbddc", "LIHTC" = "#fec44f")) +
  labs(
    title = "Median Change in Signed DI by Decade and LIHTC Treatment",
    x = "Decade",
    y = "Median ΔDI",
    fill = "LIHTC Treatment"
  ) +
  geom_text(
    aes(label = round(median_delta, 2),
        vjust = ifelse(median_delta >= 0, 1.5, -0.5)), 
    position = position_dodge(width = 0.6),
    size = 3.5,
    color = "black"
  )+
  theme_bw(base_size = 13) +  theme(legend.position = "top", legend.title = element_text(face = "bold"),
                                    plot.title = element_text(face = "bold", size = 16, margin = margin(b = 6)),
                                    plot.subtitle = element_text(size = 13, margin = margin(b = 10)),
                                    plot.caption = element_text(size=12, hjust = 0.5),
                                    axis.title.y = element_text(margin = margin(r = 10)),
                                    axis.title.x = element_text(margin = margin(t = 10))
  )

ggsave("output/Q2.2-bars.png", width = 9, height = 6, dpi = 300)
#----
#4.2.2： ever treated vs. never treated
ever_treated <- tract_wide %>%
  mutate(
    delta_80_19 = signed_percentile_2019 - signed_percentile_1980,
    ever_treated = as.integer(rowSums(select(., starts_with("lihtc_placed_")), na.rm = TRUE) > 0),
    treated = factor(ever_treated, levels = c(0,1), labels = c("Never LIHTC", "Ever LIHTC"))
  )

ggplot(ever_treated,aes(x = treated, y = delta_80_19, fill = treated)) +
  geom_boxplot(outlier.shape = NA) +
  labs(
    title = "Change in DI from 1980 to 2019 by Ever LIHTC Treatment",
    y = "Change in DI (1980–2019)", x = ""
  ) + stat_summary(
    fun = median,
    geom = "text",
    aes(label = round(..y.., 1)),
    position = position_dodge(width = 0.75),
    vjust = -0.5, size = 3.3, color = "black"
  ) +  theme_bw()+
  scale_fill_brewer(name = "LIHTC Placement Status", palette = "Set1") +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top",
        legend.title = element_text(face = "bold")) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "#ffffbf") + coord_cartesian(ylim = c(-100, 100))

ggsave("output/Q2.3-ever-never.png", width = 9, height = 6, dpi = 300)
