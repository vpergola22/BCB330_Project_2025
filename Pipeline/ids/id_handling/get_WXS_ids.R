# Load necessary library
library(dplyr)

# Read the TSV file
data <- read.csv("/ddn_exa/campbell/vpergola/Data/Tumour/wes/song2022/metadata/PRJNA754592.runinfo_ftp.tsv", sep = "\t", header = TRUE) #TODO can change this to another file in future

# Filter rows where library_strategy is "WXS" and select relevant columns
filtered_data <- data %>%
  filter(library_strategy == "WXS") %>%
  select(run_accession, experiment_accession, sample_accession)

# Rename the columns
colnames(filtered_data) <- c("SRR_num", "SRX_num", "SAMN_num")

# Save the filtered data to a CSV file
write.csv(filtered_data, "WSX_ids.csv", row.names = FALSE)

# Print a message to indicate the script has finished
cat("Filtered data saved to 'filtered_data.csv'.\n")
