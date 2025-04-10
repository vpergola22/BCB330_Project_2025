# Load required libraries
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(stringr)

# Function to process a single CSV file and count D gene presence/absence
process_csv_file <- function(file_path) {
  # Extract sample name from file path
  sample_name <- tools::file_path_sans_ext(basename(file_path))
  
  # Read the CSV file
  df <- read_csv(file_path, show_col_types = FALSE)
  
  # Check if d_gene column exists
  if (!"d_gene" %in% colnames(df)) {
    stop("Column 'd_gene' not found in file: ", file_path)
  }
  
  # Count D gene presence/absence
  d_status <- df %>%
    # Create a Has_D_Gene column based on d_gene values
    mutate(Has_D_Gene = ifelse(is.na(d_gene) | d_gene == "None" | d_gene == "", "No", "Yes")) %>%
    count(Has_D_Gene) %>%
    mutate(Sample = sample_name,
           Proportion = n / sum(n))
  
  # Make sure both "Yes" and "No" categories exist
  all_categories <- data.frame(
    Has_D_Gene = c("Yes", "No")
  )
  
  d_status <- d_status %>%
    right_join(all_categories, by = "Has_D_Gene") %>%
    mutate(
      Sample = sample_name,
      n = ifelse(is.na(n), 0, n),
      Proportion = n / sum(n, na.rm = TRUE)
    )
  
  return(d_status)
}

# Function to process multiple CSV files in a directory
process_directory <- function(directory) {
  # Get list of all CSV files
  csv_files <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(csv_files) == 0) {
    stop("No CSV files found in directory: ", directory)
  }
  
  message("Found ", length(csv_files), " CSV files")
  
  # Initialize empty dataframe to store results
  all_samples_d_status <- data.frame()
  
  # Process each file
  for (file in csv_files) {
    message("Processing file: ", basename(file))
    file_results <- process_csv_file(file)
    all_samples_d_status <- bind_rows(all_samples_d_status, file_results)
  }
  
  return(all_samples_d_status)
}

# Function to process a single CSV file with multiple samples
process_single_csv_with_barcodes <- function(file_path) {
  # Read the CSV file
  df <- read_csv(file_path, show_col_types = FALSE)
  
  # Check if required columns exist
  if (!all(c("d_gene", "barcode") %in% colnames(df))) {
    stop("Required columns 'd_gene' and 'barcode' not found in file: ", file_path)
  }
  
  # Extract sample names from barcodes (everything before the first dash)
  df <- df %>%
    mutate(Sample = str_replace(barcode, "^(.+?)-.+$", "\\1"))
  
  # Count D gene presence/absence for each sample
  d_status <- df %>%
    # Create a Has_D_Gene column based on d_gene values
    mutate(Has_D_Gene = ifelse(is.na(d_gene) | d_gene == "None" | d_gene == "", "No", "Yes")) %>%
    group_by(Sample, Has_D_Gene) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(Sample) %>%
    mutate(Proportion = n / sum(n))
  
  # Ensure both "Yes" and "No" categories exist for each sample
  samples <- unique(d_status$Sample)
  categories <- data.frame(
    Sample = rep(samples, each = 2),
    Has_D_Gene = rep(c("Yes", "No"), length(samples))
  )
  
  d_status <- categories %>%
    left_join(d_status, by = c("Sample", "Has_D_Gene")) %>%
    mutate(
      n = ifelse(is.na(n), 0, n),
      Proportion = ifelse(is.na(Proportion), 0, Proportion)
    )
  
  return(d_status)
}

#========= MAIN EXECUTION =========
# Choose which function to use based on your data structure
# If you have one CSV file per sample:
csv_dir <- "/ddn_exa/campbell/vpergola/Data/Tcell/sc/song2022/not_filtered"
all_samples_d_status <- process_directory(csv_dir)

# If you have a single CSV file with multiple samples identified by barcodes:
#csv_file <- "path/to/your/single/file.csv" # Update this path
#all_samples_d_status <- process_single_csv_with_barcodes(csv_file)

# View the results
print(all_samples_d_status)

# Plot the results - comparing Yes/No for each sample with percentages on y-axis
ggplot(all_samples_d_status, aes(x = Sample, y = Proportion, fill = Has_D_Gene)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = n), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "D Gene Presence by Sample",
       x = "Sample", 
       y = "Percentage", 
       fill = "Has D Gene")

# Plot showing only missing D genes with percentages
missing_d <- all_samples_d_status %>%
  filter(Has_D_Gene == "No") %>%
  arrange(desc(Proportion))

ggplot(missing_d, aes(x = reorder(Sample, Proportion), y = Proportion)) +
  geom_bar(stat = "identity", fill = "firebrick") +
  geom_text(aes(label = n), vjust = -0.5, size = 3) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Missing D Genes by Sample (Percentage)",
       x = "Sample", 
       y = "Percentage Missing D Genes") +
  scale_y_continuous(labels = scales::percent, 
                     limits = c(0, max(missing_d$Proportion) * 1.1)) # Room for labels