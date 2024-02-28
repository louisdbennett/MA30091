library(ggplot2)
library(GGally)
library(dplyr)

spruce_data <- read.csv('spruce1994.csv') %>% 
  mutate(
    across(c(relief, soiltxt, soiltype), stringr::str_trim),
    soiltype = if_else(grepl('brown', soiltype), 'brown soil', soiltype),
    area = case_match(
      area,
      'around Stuttgart/Neckarland' ~ 'Neckar',
      'Baar/Blackforest' ~ 'Baar',
      'Black Forest' ~ 'Black',
      'Donau, Lake of Konstanz' ~ 'Donau',
      'Oden Forest' ~ 'Oden',
      'Rhine area' ~ 'Rhine',
      'Swabian Alp' ~ 'Swabia',
      .default = area
    )
  ) %>% 
  select(-X)

unique(spruce_data$relief)

spruce_data %>% count(soiltype)

col_types <- sapply(colnames(spruce_data), function(x) mode(spruce_data[[x]]))
