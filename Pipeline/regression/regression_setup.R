# Load Libraries
library(immunarch)
library(maftools)

# ======= Formatting Data =======
# Define file paths
gsm_mapping_path <- "/ddn_exa/campbell/vpergola/Data/Pipeline/immunarch/pt_gsm_mapping.csv"
data_directory <- "/ddn_exa/campbell/vpergola/Data/Tcell/sc/song2022/immunarch"

# Step 2: get all paths to files in the directory
list_dirs <- list.files(data_directory, full.names = TRUE)

# Step 3: Load the selected files into an immdata object using repLoad
if (length(list_dirs) > 0) {
  # Load the files
  immdata <- repLoad(list_dirs)
  
  # Print a message confirming the files were loaded
  print(paste("Successfully loaded", length(list_dirs), "files into immdata."))
} else {
  print("No valid files found.")
}

#======= Get Public TCRs =======
public_tcrs <- pubRep(immdata$data, .col = "v+j+aa", .quant = "count", .coding = TRUE)

# Filter out by number of samples
public_tcrs <- public_tcrs %>% filter(Samples > 1)

public_tcrs$combined <- apply(public_tcrs, 1, function(x) paste(x, collapse = ","))

# Save to CSV
write.csv(public_tcrs, "public_tcrs_vjaa.csv", row.names = FALSE)



# ======= Load in MAF files =======
# Function to read and combine MAF files more carefully
combine_maf_files <- function(directory, clinical_file = NULL) {
  # List all .maf or .maf.gz files in the directory
  maf_files <- list.files(directory, pattern = "\\.maf(\\.gz)?$", full.names = TRUE)
  
  if(length(maf_files) == 0) {
    stop("No MAF files found in the specified directory.")
  }
  
  # Print the files found
  cat("Found", length(maf_files), "MAF files:\n")
  for(file in maf_files) {
    cat(" -", basename(file), "\n")
  }
  
  # Read all MAF files individually first to ensure they all work
  individual_mafs <- list()
  sample_ids <- character(0)
  
  for(i in seq_along(maf_files)) {
    maf_file <- maf_files[i]
    sample_name <- tools::file_path_sans_ext(basename(maf_file))
    cat("Reading file", i, "of", length(maf_files), ":", basename(maf_file), "\n")
    
    # Read the MAF file
    tryCatch({
      individual_mafs[[i]] <- read.maf(maf = maf_file)
      
      # Get sample IDs from this file and add to our tracking list
      file_samples <- getSampleSummary(individual_mafs[[i]])$Tumor_Sample_Barcode
      cat("  Sample count in this file:", length(file_samples), "\n")
      cat("  Sample IDs:", paste(file_samples, collapse=", "), "\n")
      sample_ids <- c(sample_ids, file_samples)
      
      cat("  Mutation count:", nrow(individual_mafs[[i]]@data), "\n")
      cat("  Number of columns:", ncol(individual_mafs[[i]]@data), "\n")
      
    }, error = function(e) {
      cat("ERROR reading file", basename(maf_file), ":", conditionMessage(e), "\n")
      NULL
    })
  }
  
  # Remove any NULL entries from failed reads
  individual_mafs <- individual_mafs[!sapply(individual_mafs, is.null)]
  
  # Check for duplicate sample IDs
  dup_samples <- sample_ids[duplicated(sample_ids)]
  if(length(dup_samples) > 0) {
    cat("WARNING: Duplicate sample IDs found:", paste(dup_samples, collapse=", "), "\n")
    cat("This could cause samples to be overwritten during merging.\n")
  }
  
  # Method 1: Merge at the raw data level
  cat("\nAttempting to merge MAF files at the data level...\n")
  
  # Get raw data from each MAF object and ensure consistent columns
  all_columns <- unique(unlist(lapply(individual_mafs, function(maf) colnames(maf@data))))
  cat("Total unique columns across all MAF files:", length(all_columns), "\n")
  
  # Prepare a list to store the standardized data frames
  standardized_data <- list()
  
  # Standardize each data frame to have all columns
  for(i in seq_along(individual_mafs)) {
    df <- individual_mafs[[i]]@data
    
    # Find missing columns
    missing_cols <- setdiff(all_columns, colnames(df))
    
    if(length(missing_cols) > 0) {
      cat("  File", i, "is missing", length(missing_cols), "columns.\n")
      
      # Add missing columns with NA values
      for(col in missing_cols) {
        df[[col]] <- NA
      }
    }
    
    # Store the standardized data frame (using data.table syntax)
    standardized_data[[i]] <- df[, ..all_columns]
    cat("  Processed file", i, "with dimensions:", dim(standardized_data[[i]]), "\n")
  }
  
  # Combine all data frames
  combined_data <- do.call(rbind, standardized_data)
  cat("Combined data dimensions:", dim(combined_data), "\n")
  
  # Read clinical data if provided
  clinical_data <- NULL
  if(!is.null(clinical_file) && file.exists(clinical_file)) {
    cat("Reading clinical data from:", clinical_file, "\n")
    clinical_data <- read.delim(clinical_file, stringsAsFactors = FALSE)
  }
  
  # Create the final MAF object
  cat("Creating final MAF object...\n")
  combined_maf <- read.maf(maf = combined_data, clinicalData = clinical_data)
  
  # Verify all samples are present
  samples_in_result <- getSampleSummary(combined_maf)$Tumor_Sample_Barcode
  cat("Samples in final MAF object:", length(samples_in_result), "\n")
  cat("Sample IDs:", paste(samples_in_result, collapse=", "), "\n")
  
  # Check if any samples were lost
  missing_samples <- setdiff(sample_ids, samples_in_result)
  
  if(length(missing_samples) > 0) {
    cat("WARNING: The following samples are missing from the final result:", 
        paste(missing_samples, collapse=", "), "\n")
  }
  
  return(combined_maf)
}

# Set the directory path
maf_directory <- "/ddn_exa/campbell/vpergola/Data/Pipeline/meskit/final_mafs"

# Set the clinical data file path
clinical_file <- file.path(maf_directory, "clinical_data.tsv")

# Process MAF files
cat("Starting MAF file processing...\n")
combined_maf <- combine_maf_files(maf_directory, clinical_file)


#======= Tumour Mutation Frequencies =======
mutation_counts <- table(combined_maf@data$Hugo_Symbol)

