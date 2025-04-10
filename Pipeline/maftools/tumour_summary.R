#!/usr/bin/R

# ======= Libraries =======
library(maftools)
library(dplyr)


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

# Basic exploration of the MAF object
print(combined_maf)

# Get sample summary
sample_summary <- getSampleSummary(combined_maf)
print("Sample Summary:")
print(sample_summary)

# Get gene summary
gene_summary <- getGeneSummary(combined_maf)
print("Gene Summary:")
print(head(gene_summary))

# Optional: Write MAF summary to file
write.mafSummary(maf = combined_maf, basename = "combined_maf_summary")


# ======= Analysis =======
# TODO FIXES
# ======= Sample Labeling Setup =======
cat("Setting up sample labeling from clinical data...\n")

# Ensure clinical data is properly loaded and matched
if(!is.null(clinical_file) && file.exists(clinical_file)) {
  # Read the clinical data
  clinical_data <- read.delim(clinical_file, stringsAsFactors = FALSE)
  
  # Check if 'patient_ID' column exists in clinical data
  if("Patient_ID" %in% colnames(clinical_data)) {
    cat("Found Patient_ID column in clinical data\n")
    
    # Get current sample IDs from MAF
    current_samples <- getSampleSummary(combined_maf)$Tumor_Sample_Barcode
    
    # Create a mapping table between sample IDs and patient IDs
    # Assuming the column that matches sample IDs in clinical data is 'Tumor_Sample_Barcode'
    # If it's different in your clinical data, change it accordingly
    sample_col <- if("Tumor_Sample_Barcode" %in% colnames(clinical_data)) 
                    "Tumor_Sample_Barcode" else colnames(clinical_data)[1]
    
    cat("Using", sample_col, "from clinical data to match with MAF sample IDs\n")
    
    # Create the mapping
    sample_to_patient <- clinical_data[, c(sample_col, "Patient_ID")]
    colnames(sample_to_patient)[1] <- "Tumor_Sample_Barcode"
    
    # Check if all MAF samples are in the clinical data
    missing_samples <- setdiff(current_samples, sample_to_patient$Tumor_Sample_Barcode)
    if(length(missing_samples) > 0) {
      cat("WARNING: Some samples in MAF are not in clinical data:", 
          paste(missing_samples, collapse=", "), "\n")
    }
    
    # Create a function to generate custom sample labels
    get_patient_labels <- function(samples) {
      # Match samples to patient IDs
      patient_ids <- sample_to_patient$Patient_ID[match(samples, sample_to_patient$Tumor_Sample_Barcode)]
      # Replace NA with original sample ID
      patient_ids[is.na(patient_ids)] <- samples[is.na(patient_ids)]
      return(patient_ids)
    }
    
    cat("Sample to patient ID mapping set up successfully\n")
  } else {
    cat("WARNING: 'Patient_ID' column not found in clinical data. Will use original sample IDs.\n")
    # Create a dummy function that returns the original sample IDs
    get_patient_labels <- function(samples) { return(samples) }
  }
} else {
  cat("WARNING: Clinical data file not found or not specified. Will use original sample IDs.\n")
  # Create a dummy function that returns the original sample IDs
  get_patient_labels <- function(samples) { return(samples) }
}

# ======= Create Patient-Labeled MAF Object =======
# This function creates a new MAF object with patient IDs instead of sample IDs
create_patient_labeled_maf <- function(maf_obj, clinical_data) {
  # Check if we have the necessary columns
  if(!("Tumor_Sample_Barcode" %in% colnames(clinical_data)) || 
     !("Patient_ID" %in% colnames(clinical_data))) {
    cat("WARNING: Cannot create patient-labeled MAF: missing required columns in clinical data\n")
    return(maf_obj)  # Return original if we can't proceed
  }
  
  # Create a mapping dictionary
  sample_to_patient <- setNames(
    clinical_data$Patient_ID, 
    clinical_data$Tumor_Sample_Barcode
  )
  
  # Make a copy of the MAF data
  new_maf_data <- maf_obj@data
  
  # Replace sample IDs with patient IDs
  new_maf_data$Tumor_Sample_Barcode <- sapply(
    new_maf_data$Tumor_Sample_Barcode,
    function(s) {
      if(s %in% names(sample_to_patient)) {
        return(sample_to_patient[s])
      } else {
        return(s)  # Keep original if no mapping exists
      }
    }
  )
  
  # Create a new MAF object with the modified data
  patient_maf <- read.maf(maf = new_maf_data, clinicalData = maf_obj@clinical.data)
  
  return(patient_maf)
}

# Create a patient-labeled MAF object if possible
patient_labeled_maf <- NULL
if(!is.null(clinical_file) && file.exists(clinical_file) && 
   "Patient_ID" %in% colnames(clinical_data) && 
   "Tumor_Sample_Barcode" %in% colnames(clinical_data)) {
  
  cat("Creating a MAF object with patient IDs instead of sample IDs...\n")
  patient_labeled_maf <- create_patient_labeled_maf(combined_maf, clinical_data)
}

# Get the top mutated samples for rainfall plots
sample_summary <- getSampleSummary(combined_maf)
# Sort by mutation count in descending order
sample_summary <- sample_summary[order(-sample_summary$total), ]
# Take top 5 samples (or adjust as needed)
top_mutated_samples <- sample_summary$Tumor_Sample_Barcode[1:min(5, nrow(sample_summary))]

cat("Top mutated samples for rainfall plots:", paste(top_mutated_samples, collapse=", "), "\n")

# Define top genes for lollipop plots
gene_summary <- getGeneSummary(combined_maf)
top_genes <- gene_summary$Hugo_Symbol[1:min(5, nrow(gene_summary))]

# ======= Modified Visualization Code =======

# 1. Oncoplots with patient IDs
cat("Creating oncoplots with patient IDs...\n")

# Use patient-labeled MAF if available, otherwise use original MAF
maf_for_plots <- if(!is.null(patient_labeled_maf)) patient_labeled_maf else combined_maf

# Standard oncoplot
pdf("oncoplot_patient_labeled.pdf", width=12, height=10)
# Check if we have data before plotting
if(nrow(maf_for_plots@data) > 0) {
  oncoplot(maf = maf_for_plots, 
         top = 20,
         fontSize = 0.8,
         sampleOrder = getSampleSummary(maf_for_plots)$Tumor_Sample_Barcode,
         showTumorSampleBarcodes = TRUE,
         barcode_mar = 6,  
         gene_mar = 5)
} else {
  # Create a blank plot with a message if no data
  plot(1, type="n", axes=FALSE, xlab="", ylab="")
  text(1, 1, "No mutation data available for oncoplot", cex=1.5)
}
dev.off()

# For GISTIC data (if available)
if(exists("gistic_attempt") && !is.null(gistic_attempt)) {
  pdf("oncoplot_CN_patient_labeled.pdf", width=14, height=10)
  if(nrow(maf_for_plots@data) > 0) {
    oncoplot(maf = maf_for_plots,
           top = 20,
           gistic = gistic_attempt,
           fontSize = 0.8,
           sampleOrder = getSampleSummary(maf_for_plots)$Tumor_Sample_Barcode,
           showTumorSampleBarcodes = TRUE,
           barcode_mar = 6,
           gene_mar = 5)
  } else {
    plot(1, type="n", axes=FALSE, xlab="", ylab="")
    text(1, 1, "No mutation data available for GISTIC oncoplot", cex=1.5)
  }
  dev.off()
}

# Optional: Visualizations
# Oncoplot of top 10 mutated genes
pdf("oncoplot_top10.pdf", width=10, height=8)
oncoplot(maf = combined_maf, 
        top = 10,
        showTumorSampleBarcodes = TRUE,
        )
dev.off()

# Variant classification plot
pdf("variant_classification_summary.pdf", width=10, height=6)
plotmafSummary(maf = combined_maf, rmOutlier = TRUE, addStat = 'median', showBarcodes = TRUE)
dev.off()

# 2. Rainfall plots with patient IDs in filenames and titles
cat("Creating rainfall plots for top mutated samples...\n")

# Check if top_mutated_samples has any samples
if(length(top_mutated_samples) > 0) {
  for(i in 1:length(top_mutated_samples)) {
    sample <- top_mutated_samples[i]
    patient_id <- get_patient_labels(sample)
    
    output_file <- paste0("rainfall_plot_", gsub("[^a-zA-Z0-9]", "_", patient_id), ".pdf")
    cat("  Creating rainfall plot for", patient_id, "...\n")
    
    # Check if this sample has mutations
    sample_muts <- combined_maf@data[combined_maf@data$Tumor_Sample_Barcode == sample, ]
    
    if(nrow(sample_muts) > 0) {
      tryCatch({
        pdf(output_file, width=12, height=6)
        # Add a title with the patient ID after the plot
        par(mar=c(5, 4, 4, 2) + 0.1)  # Increase top margin for title
        rainfallPlot(maf = combined_maf, 
                    detectChangePoints = TRUE, 
                    pointSize = 0.4,
                    tsb = sample)
        # Add title manually after plot is created
        mtext(paste("Rainfall Plot for Patient:", patient_id), side=3, line=1, cex=1.2)
        dev.off()
        cat("    Plot created successfully.\n")
      }, error = function(e) {
        cat("    Error creating rainfall plot:", conditionMessage(e), "\n")
        # Make sure PDF device is closed in case of error
        if(dev.cur() > 1) dev.off()
      })
    } else {
      cat("    No mutations found for this sample. Skipping plot.\n")
    }
  }
} else {
  cat("No samples found for rainfall plots. Skipping this section.\n")
}

# 3. VAF Plot with patient IDs
# Check if the vaf_col_to_use variable exists
vaf_col_to_use <- if(exists("vaf_col_to_use")) vaf_col_to_use else "i_tumor_f"

# Check if the column exists in our data
if(vaf_col_to_use %in% colnames(combined_maf@data)) {
  pdf("vaf_boxplot.pdf", width=12, height=8)
  # Get the original sample order
  original_samples <- unique(combined_maf@data$Tumor_Sample_Barcode)
  
  # Plot without custom labels (not supported in this version)
  plotVaf(maf = combined_maf, 
          vafCol = vaf_col_to_use,
          sampleOrder = original_samples)
  dev.off()
}

# 4. Create a custom summary table with patient IDs
cat("Creating patient ID summary table...\n")
sample_summary <- getSampleSummary(combined_maf)

# Create a standard data frame (not a data.table) to avoid issues
sample_summary_df <- as.data.frame(sample_summary)
# Add the patient ID column
sample_summary_df$patient_ID <- get_patient_labels(sample_summary_df$Tumor_Sample_Barcode)
write.table(sample_summary_df, "patient_sample_summary.tsv", sep="\t", quote=FALSE, row.names=FALSE)


# 6. Somatic Interactions
cat("Creating somatic interactions plot...\n")
tryCatch({
  pdf("somatic_interactions.pdf", width=10, height=10)
  # No title parameter in this version
  somaticInteractions(maf = combined_maf, 
                     top = 25,
                     pvalue = 0.05)
  # Add title manually
  mtext("Somatic Interactions in Patient Cohort", side=3, line=1, cex=1.2)
  dev.off()
}, error = function(e) {
  cat("Error creating somatic interactions plot:", conditionMessage(e), "\n")
  if(dev.cur() > 1) dev.off()
})

# 7. Oncogenic Pathways
tryCatch({
  pdf("oncogenic_pathways.pdf", width=12, height=16)
  # No title parameter in this version
  pathways(maf = combined_maf)
  # Add title manually
  mtext("Oncogenic Pathways in Patient Cohort", side=3, line=1, cex=1.2)
  dev.off()
}, error = function(e) {
  cat("Error creating oncogenic pathways plot:", conditionMessage(e), "\n")
  if(dev.cur() > 1) dev.off()
})

# 8. Lollipop plots with patient mutation details
# Determine a valid AACol for lollipop plots
aa_col_to_use <- if(exists("aa_col_to_use")) aa_col_to_use else "Protein_Change"

if(aa_col_to_use %in% colnames(combined_maf@data) && length(top_genes) > 0) {
  for(gene in top_genes) {
    output_file <- paste0("lollipop_", gene, ".pdf")
    
    # Try to create the lollipop plot
    tryCatch({
      pdf(output_file, width=12, height=6)
      # Check if lollipopPlot supports repel parameter
      has_repel <- "repel" %in% names(formals(lollipopPlot))
      
      # Call with appropriate parameters for this version
      if(has_repel) {
        lollipopPlot(maf = combined_maf, 
                    gene = gene, 
                    AACol = aa_col_to_use,
                    showMutationRate = TRUE,
                    labelPos = NULL,
                    repel = TRUE)
      } else {
        lollipopPlot(maf = combined_maf, 
                    gene = gene, 
                    AACol = aa_col_to_use,
                    showMutationRate = TRUE)
      }
      
      # Add title manually
      mtext(paste(gene, "Mutations in Patient Cohort"), side=3, line=1, cex=1.2)
      dev.off()
      cat("Created lollipop plot for", gene, "\n")
    }, error = function(e) {
      cat("Error creating lollipop plot for", gene, ":", conditionMessage(e), "\n")
      if(dev.cur() > 1) dev.off()
    })
  }
}

# ======= Additional Analysis for Variable Genes and TCR Presence =======

# Create oncoplot for genes mutated in at least 50% of samples
cat("Creating oncoplot for genes mutated in at least 50% of samples...\n")

# Calculate the number of samples
total_samples <- length(unique(getSampleSummary(combined_maf)$Tumor_Sample_Barcode))
threshold_50_percent <- ceiling(total_samples * 0.5)

# Get genes mutated in at least 50% of samples
gene_summary <- getGeneSummary(combined_maf)
genes_50_percent <- gene_summary$Hugo_Symbol[gene_summary$MutatedSamples >= threshold_50_percent]

# If there are genes to plot, split them into chunks of 100 and create multiple plots
if(length(genes_50_percent) > 0) {
  # Calculate how many plots we need
  num_plots <- ceiling(length(genes_50_percent) / 100)
  cat(paste0("Creating ", num_plots, " oncoplots with up to 100 genes each...\n"))
  
  # Split genes into chunks of 100
  for(i in 1:num_plots) {
    # Calculate start and end indices for this chunk
    start_idx <- ((i-1) * 100) + 1
    end_idx <- min(i * 100, length(genes_50_percent))
    
    # Get the genes for this chunk
    current_genes <- genes_50_percent[start_idx:end_idx]
    
    # Create a PDF for this chunk
    pdf_filename <- paste0("oncoplot_50percent_part", i, "_of_", num_plots, ".pdf")
    pdf(pdf_filename, width=12, height=10)
    
    cat(paste0("Creating plot ", i, " with genes ", start_idx, " to ", end_idx, "...\n"))
    
    oncoplot(maf = maf_for_plots, 
             genes = current_genes,
             fontSize = 0.8,
             showTumorSampleBarcodes = TRUE,
             barcode_mar = 6,
             gene_mar = 5,
             titleText = paste0("Oncoplot Part ", i, " of ", num_plots, " (Genes mutated in ≥50% of samples)")
             )  # Include clinical annotations if available
    
    dev.off()
    cat(paste0("Saved ", pdf_filename, "\n"))
  }
} else {
  # Create a single plot indicating no genes found
  pdf("oncoplot_50percent_mutated.pdf", width=12, height=10)
  plot(1, type="n", axes=FALSE, xlab="", ylab="")
  text(1, 1, "No genes found mutated in at least 50% of samples", cex=1.5)
  dev.off()
  cat("No genes found mutated in at least 50% of samples\n")
}

# TODO do multiple plots with max 100 genes (in a loop) 



# ======= Analyzing Mutation Distribution Statistics =======

# Generate a detailed report of mutation percentages
cat("Generating mutation percentage report...\n")

# Function to calculate and explain mutation percentages
analyze_mutation_percentages <- function(maf_obj) {
  # Get gene summary
  gene_summ <- getGeneSummary(maf_obj)
  
  # Get sample count
  total_samples <- length(unique(getSampleSummary(maf_obj)$Tumor_Sample_Barcode))
  
  # Calculate percentage of samples with mutations for each gene
  gene_summ$MutationPercentage <- (gene_summ$MutatedSamples / total_samples) * 100
  
  # Sort by mutation percentage (descending)
  gene_summ <- gene_summ[order(-gene_summ$MutationPercentage), ]
  
  # Create groups based on mutation percentage
  gene_summ$MutationGroup <- cut(gene_summ$MutationPercentage, 
                                breaks = c(0, 25, 50, 75, 100), 
                                labels = c("0-25%", "26-50%", "51-75%", "76-100%"),
                                include.lowest = TRUE)
  
  # Count genes in each group
  group_counts <- table(gene_summ$MutationGroup)
  
  # Write the detailed report
  write.csv(gene_summ, "gene_mutation_percentages.csv", row.names = FALSE)
  
  # Create a summary report
  sink("mutation_percentage_summary.txt")
  cat("MUTATION PERCENTAGE ANALYSIS\n")
  cat("============================\n\n")
  cat("Total samples analyzed:", total_samples, "\n\n")
  cat("Genes by mutation frequency:\n")
  cat("- 100% of samples:", sum(gene_summ$MutationPercentage == 100), "genes\n")
  cat("- 75-99% of samples:", sum(gene_summ$MutationPercentage >= 75 & gene_summ$MutationPercentage < 100), "genes\n")
  cat("- 50-74% of samples:", sum(gene_summ$MutationPercentage >= 50 & gene_summ$MutationPercentage < 75), "genes\n")
  cat("- 25-49% of samples:", sum(gene_summ$MutationPercentage >= 25 & gene_summ$MutationPercentage < 50), "genes\n")
  cat("- 1-24% of samples:", sum(gene_summ$MutationPercentage > 0 & gene_summ$MutationPercentage < 25), "genes\n\n")
  
  cat("Top 20 most frequently mutated genes:\n")
  top20 <- head(gene_summ, 20)
  for(i in 1:nrow(top20)) {
    cat(sprintf("%d. %s: %.1f%% (%d/%d samples)\n", 
                i, top20$Hugo_Symbol[i], top20$MutationPercentage[i], 
                top20$MutatedSamples[i], total_samples))
  }
  sink()
  
  return(gene_summ)
}

# Run the analysis
mutation_percentage_analysis <- analyze_mutation_percentages(combined_maf)

# ======= TCR Presence Analysis =======

# Function to identify and analyze potential TCR-related genes
analyze_tcr_genes <- function(maf_obj) {
  # List of TCR-related genes to check
  tcr_genes <- c(
    # TCR alpha/delta genes
    "TRAC", "TRAD", "TRAJ", "TRAV",
    # TCR beta genes
    "TRBC", "TRBD", "TRBJ", "TRBV",
    # TCR gamma genes
    "TRGC", "TRGJ", "TRGV",
    # Additional TCR-related genes
    "CD3D", "CD3E", "CD3G", "CD247",
    # Check for any gene starting with TR
    grep("^TR[ABGD][CVJ]", getGeneSummary(maf_obj)$Hugo_Symbol, value = TRUE)
  )
  
  # Get unique TCR genes that exist in our data
  gene_summ <- getGeneSummary(maf_obj)
  tcr_genes_present <- intersect(unique(tcr_genes), gene_summ$Hugo_Symbol)
  
  if(length(tcr_genes_present) > 0) {
    # Create an oncoplot specifically for TCR genes
    pdf("tcr_genes_oncoplot.pdf", width=12, height=8)
    
    # Using only well-documented parameters for oncoplot
    oncoplot(maf = maf_obj, 
             genes = tcr_genes_present,
             fontSize = 0.8,
             showTumorSampleBarcodes = TRUE,
             titleFontSize = 1.2,
             title = "TCR-related Genes Mutation Profile")
    
    dev.off()
    
    # Create lollipop plots for each TCR gene
    if(length(tcr_genes_present) > 0) {
      pdf("tcr_genes_lollipop_plots.pdf", width=10, height=8)
      
      for(gene in tcr_genes_present) {
        # Check if there are enough mutations for this gene to create a plot
        gene_muts <- length(which(maf_obj@data$Hugo_Symbol == gene))
        
        if(gene_muts >= 3) {  # Only create lollipop plot if there are at least 3 mutations
          tryCatch({
            lollipopPlot(maf = maf_obj, 
                         gene = gene, 
                         showMutationRate = TRUE,
                         labelPos = "all",  # Label all positions
                         repel = TRUE,      # Avoid label overlapping
                         domainLabelSize = 3,
                         pointSize = 1.2,
                         axisTextSize = 0.8)
          }, error = function(e) {
            message(paste("Could not create lollipop plot for gene", gene, ":", e$message))
          })
        } else {
          message(paste("Skipping lollipop plot for", gene, "- insufficient mutations (", gene_muts, ")"))
        }
      }
      
      dev.off()
    }
    
    # Extract TCR gene data for detailed analysis
    tcr_data <- gene_summ[gene_summ$Hugo_Symbol %in% tcr_genes_present, ]
    write.csv(tcr_data, "tcr_genes_analysis.csv", row.names = FALSE)
    
    # Create a summary report
    sink("tcr_analysis_summary.txt")
    cat("TCR GENE MUTATION ANALYSIS\n")
    cat("==========================\n\n")
    cat("TCR-related genes identified in dataset:", length(tcr_genes_present), "\n\n")
    
    # List all TCR genes found
    cat("TCR genes present in dataset:\n")
    for(i in 1:length(tcr_genes_present)) {
      gene_row <- tcr_data[tcr_data$Hugo_Symbol == tcr_genes_present[i], ]
      cat(sprintf("%d. %s: Mutated in %d samples (%.1f%%)\n", 
                  i, tcr_genes_present[i], gene_row$MutatedSamples, 
                  (gene_row$MutatedSamples/length(unique(getSampleSummary(maf_obj)$Tumor_Sample_Barcode)))*100))
    }
    sink()
    
    return(tcr_data)
  } else {
    cat("No TCR-related genes found in the dataset.\n")
    return(NULL)
  }
}

# Run TCR analysis
tcr_analysis <- analyze_tcr_genes(combined_maf)

# ======= Additional Visualization for Understanding 100% Mutations =======

# Create a heatmap showing mutation types for 100% mutated genes
cat("Creating heatmap of mutation types for genes mutated in all samples...\n")

# Get genes mutated in 100% of samples
total_samples <- length(unique(getSampleSummary(combined_maf)$Tumor_Sample_Barcode))
genes_100percent <- mutation_percentage_analysis$Hugo_Symbol[
  mutation_percentage_analysis$MutatedSamples == total_samples
]

# Generate mutation type summary for these genes
if(length(genes_100percent) > 0) {
  # Create a detailed mutation type breakdown
  mutation_type_summary <- NULL
  for(gene in genes_100percent) {
    gene_muts <- subsetMaf(combined_maf, genes = gene)
    if(!is.null(gene_muts)) {
      mut_types <- table(gene_muts@data$Variant_Classification)
      mut_types_df <- data.frame(
        Gene = gene,
        Mutation_Type = names(mut_types),
        Count = as.numeric(mut_types)
      )
      mutation_type_summary <- rbind(mutation_type_summary, mut_types_df)
    }
  }
  
  # Save this data
  if(!is.null(mutation_type_summary)) {
    write.csv(mutation_type_summary, "mutation_types_100percent_genes.csv", row.names = FALSE)
  }
  
  # Create proper visualization with a maximum of 100 genes per plot
  if(length(genes_100percent) > 0) {
    # Calculate how many plots we need
    max_genes_per_plot <- 100
    num_plots <- ceiling(length(genes_100percent) / max_genes_per_plot)
    cat(paste0("Creating ", num_plots, " oncoplots with up to 100 genes each...\n"))
    
    # Split genes into chunks of 100
    for(i in 1:num_plots) {
      # Calculate start and end indices for this chunk
      start_idx <- ((i-1) * max_genes_per_plot) + 1
      end_idx <- min(i * max_genes_per_plot, length(genes_100percent))
      
      # Get the genes for this chunk
      current_genes <- genes_100percent[start_idx:end_idx]
      
      # Create a PDF for this chunk
      pdf_filename <- paste0("mutation_types_100percent_genes_part", i, "_of_", num_plots, ".pdf")
      pdf(pdf_filename, width=12, height=10)
      
      cat(paste0("Creating plot ", i, " with genes ", start_idx, " to ", end_idx, "...\n"))
      
      # Use oncoplot for this subset of genes
      oncoplot(maf = combined_maf, 
               genes = current_genes,
               showTumorSampleBarcodes = TRUE,
               removeNonMutated = FALSE,
               titleText = paste0("Oncoplot Part ", i, " of ", num_plots, " (Genes mutated in 100% of samples)"))
      
      dev.off()
      cat(paste0("Saved ", pdf_filename, "\n"))
    }
    
#     # Create lollipop plots in a separate file if there are not too many genes
#     if(length(genes_100percent) <= 10) {
#       pdf_filename <- "lollipop_plots_100percent_genes.pdf"
#       pdf(pdf_filename, width=12, height=10)
#       for(gene in genes_100percent) {
#         lollipopPlot(maf = combined_maf, gene = gene, showMutationRate = TRUE)
#       }
#       dev.off()
#       cat(paste0("Saved lollipop plots to ", pdf_filename, "\n"))
#     } else {
#       # Create multiple files for lollipop plots if needed
#       lollipop_batch_size <- 10
#       num_lollipop_files <- ceiling(length(genes_100percent) / lollipop_batch_size)
      
#       for(j in 1:num_lollipop_files) {
#         start_idx <- ((j-1) * lollipop_batch_size) + 1
#         end_idx <- min(j * lollipop_batch_size, length(genes_100percent))
#         current_lollipop_genes <- genes_100percent[start_idx:end_idx]
        
#         pdf_filename <- paste0("lollipop_plots_100percent_genes_part", j, "_of_", num_lollipop_files, ".pdf")
#         pdf(pdf_filename, width=12, height=10)
        
#         for(gene in current_lollipop_genes) {
#           lollipopPlot(maf = combined_maf, gene = gene, showMutationRate = TRUE)
#         }
        
#         dev.off()
#         cat(paste0("Saved lollipop plots to ", pdf_filename, "\n"))
#       }
#     }
   }
 }

# ======= Top Mutated Genes Lollipop Plot Analysis =======

# Function to create lollipop plots for top mutated genes
analyze_top_mutated_genes <- function(maf_obj, n_genes = 15, min_mutations = 5) {
  # Get the summary of genes and sort by mutation frequency
  gene_summary <- getGeneSummary(maf_obj)
  gene_summary <- gene_summary[order(gene_summary$MutatedSamples, decreasing = TRUE), ]
  
  # Select top N genes with at least min_mutations
  top_genes <- gene_summary$Hugo_Symbol[gene_summary$MutatedSamples >= min_mutations]
  if(length(top_genes) > n_genes) {
    top_genes <- top_genes[1:n_genes]
  }
  
  # Initialize top_genes_summary variable to avoid error if all plots fail
  top_genes_summary <- gene_summary[gene_summary$Hugo_Symbol %in% top_genes, ]
  
  if(length(top_genes) > 0) {
    cat(paste("\nCreating lollipop plots for top", length(top_genes), "mutated genes...\n"))
    
    # Create a PDF for all lollipop plots
    pdf("top_mutated_genes_lollipop_plots.pdf", width=10, height=8)
    
    # Create a counter for successful plots
    successful_plots <- 0
    
    for(gene in top_genes) {
      # Print gene info for tracking
      gene_row <- gene_summary[gene_summary$Hugo_Symbol == gene, ]
      cat(sprintf("Processing %s: %d mutations in %d samples\n", 
                  gene, gene_row$Mutations, gene_row$MutatedSamples))
      
      # Try to create the lollipop plot - FIX: remove problematic axisTextSize parameter
      tryCatch({
        lollipopPlot(maf = maf_obj, 
                   gene = gene, 
                   showMutationRate = TRUE,
                   labelPos = "all",  
                   repel = TRUE,      
                   domainLabelSize = 3,
                   pointSize = 1.2,
                   title = paste0(gene, " Mutation Profile"))
        
        successful_plots <- successful_plots + 1
      }, error = function(e) {
        cat(paste("  Error creating lollipop plot for", gene, ":", e$message, "\n"))
      })
    }
    
    dev.off()
    
    # If no plots were successful, remove the empty PDF
    if(successful_plots == 0) {
      file.remove("top_mutated_genes_lollipop_plots.pdf")
      cat("No successful lollipop plots were created. Check your mutation data.\n")
    } else {
      cat(paste("Successfully created", successful_plots, "lollipop plots out of", 
                length(top_genes), "attempted.\n"))
      
      # Create a summary table of the top mutated genes
      write.csv(top_genes_summary, "top_mutated_genes_summary.csv", row.names = FALSE)
      
      # Create oncoplot for top genes
      pdf("top_mutated_genes_oncoplot.pdf", width=12, height=8)
      oncoplot(maf = maf_obj, 
               genes = top_genes,
               fontSize = 0.8,
               showTumorSampleBarcodes = TRUE,
               title = "Top Mutated Genes Mutation Profile")
      dev.off()
    }
    
    return(top_genes_summary)
  } else {
    cat("No genes with sufficient mutations found for lollipop plots.\n")
    return(NULL)
  }
}

# Also create MAF summary and plots
create_maf_summary <- function(maf_obj) {
  # Create overall MAF summary
  pdf("maf_summary_plots.pdf", width=12, height=10)
  
  # Overall summary plot
  plotmafSummary(maf = maf_obj, addStat = 'median', dashboard = TRUE, 
                 titvRaw = TRUE, rmOutlier = TRUE)
  
  # Transition/Transversion plot
  titv = titv(maf = maf_obj, plot = FALSE)
  plotTiTv(res = titv)
  
  # Mutation type summary
  oncoplot(maf = maf_obj, 
           top = 20,
           fontSize = 0.8,
           showTumorSampleBarcodes = TRUE,
           title = "Top 20 Mutated Genes")
  
  # Try mutual exclusivity analysis on top genes
  tryCatch({
    somaticInteractions(maf = maf_obj, 
                        top = 25, 
                        pvalue = 0.05)
  }, error = function(e) {
    cat("Could not create somatic interactions plot:", e$message, "\n")
  })
  
  dev.off()
}

# Execute the analysis for top mutated genes
cat("\n\n==== RUNNING ADDITIONAL MUTATION ANALYSIS ====\n")

# Try running with error handling
tryCatch({
  top_genes_data <- analyze_top_mutated_genes(combined_maf, n_genes = 15, min_mutations = 5)
  create_maf_summary(combined_maf)
  cat("\nAnalysis complete. Check the output files for results.\n")
}, error = function(e) {
  cat("\nError in analysis:", e$message, "\n")
  cat("Trying alternative approach with minimal parameters...\n")
  
  # Try creating simplified lollipop plots with minimal parameters
  pdf("simplified_lollipop_plots.pdf", width=10, height=8)
  
  gene_summary <- getGeneSummary(combined_maf)
  top_genes <- gene_summary$Hugo_Symbol[order(gene_summary$MutatedSamples, decreasing = TRUE)][1:10]
  
  for(gene in top_genes) {
    tryCatch({
      lollipopPlot(maf = combined_maf, gene = gene)
    }, error = function(e) {
      cat(paste("Could not create simplified plot for", gene, "\n"))
    })
  }
  
  dev.off()
})

cat("All analyses completed!\n")
