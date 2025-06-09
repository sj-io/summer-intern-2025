library(tidyverse)
library(tidycensus)
census_api_key("a71f148809127ef86df2c9641d9dd7d8bb995d11", install = TRUE)

v23 <- load_variables(2023, "acs1")

tenure_21 <- get_acs("state", table = "B25003", year = 2022, survey = "acs1") |> left_join(v23, by = c("variable" = "name"))

View(tenure_21)
