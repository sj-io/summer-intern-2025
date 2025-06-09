library(tidyverse)
library(tidycensus)

v23 <- load_variables(2023, "acs1")

tenure_21 <- get_acs("county", table = "B25003", year = 2023, survey = "acs1") |> left_join(v23, by = c("variable" = "name"))
