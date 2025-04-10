library(dplyr)
library(tidyverse)

data1 <- read.csv("/ddn_exa/campbell/vpergola/Data/Tumour/wes/song2022/metadata/WSX_ids.csv")
data2 <- read.csv("/ddn_exa/campbell/vpergola/Data/Tumour/wes/song2022/metadata/GSM_ids.csv")

joined_data <- data1 %>%
  inner_join(data2, by = "SAMN_num")

vdj_set <- read.csv("/ddn_exa/campbell/vpergola/Data/Tests/song2022_vdj.csv") %>% select(GSM_num)

# Perform inner join to keep only matching GSM_num values
filtered_ids <- joined_data %>%
  semi_join(vdj_set, by = "GSM_num")  # Keeps only rows where GSM_num is in vdj

# Save the result to a new CSV file
write.csv(filtered_ids, "filtered_ids.csv", row.names = FALSE)