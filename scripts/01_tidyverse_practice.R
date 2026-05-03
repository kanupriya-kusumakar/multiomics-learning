install.packages(c("tidyverse", "palmerpenguins"))
library(tidyverse)
library(palmerpenguins)
glimpse(penguins)

penguins |>
  filter(!is.na(bill_length_mm)) |>
  group_by(species, island) |>
  summarise(mean_bill = mean(bill_length_mm), n = n()) 

# Suppress grouping message cleanly - 
penguins |>
  filter(!is.na(bill_length_mm)) |>
  group_by(species, island) |>
  summarise(mean_bill = mean(bill_length_mm), n = n(), .groups = "drop")

# Mutate & Rank
penguins |>
  filter(!is.na(bill_length_mm), !is.na(bill_depth_mm)) |>
  mutate(bill_ratio = bill_length_mm / bill_depth_mm) |>
  group_by(species) |>
  summarise(median_ratio = median(bill_ratio), .groups = "drop") |>
  arrange(desc(median_ratio))

# Pivot to long format
penguins_long <- penguins |>
  filter(!is.na(species)) |>
  pivot_longer(
    cols = where(is.numeric),
                 names_to = "measurement",
                 values_to = "value" )
glimpse(penguins_long) 

penguins_long <- penguins |>
  filter(!is.na(species)) |>
  pivot_longer(
    cols = c(bill_length_mm, bill_depth_mm, flipper_length_mm, body_mass_g),
    names_to = "measurement",
    values_to = "value"
  )
  
# Reusable function 
summarise_by_group <- function(df, group_col, value_col) {
  df |>
    filter(!is.na({{value_col}})) |>
    group_by({{group_col}}) |>
    summarise(
      mean = mean({{value_col}}),
      sd = sd({{value_col}}),
      n = n(),
      .groups = "drop"
    )
}

summarise_by_group(penguins, species, bill_length_mm)
summarise_by_group(penguins, island, body_mass_g)
