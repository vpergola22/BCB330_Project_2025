library(tidyverse)

# Load the IDs dataset
ids <- read_csv("ids.csv") #TODO edit this so that it is generalized

# Load the VDJ dataset (only GSM_num column is needed)
vdj <- read_csv("vdj.csv") %>% select(GSM_num)

# Perform inner join to keep only matching GSM_num values
filtered_ids <- ids %>%
  semi_join(vdj, by = "GSM_num")  # Keeps only rows where GSM_num is in vdj

# Save the filtered dataset
write_csv(filtered_ids, "filtered_ids.csv")