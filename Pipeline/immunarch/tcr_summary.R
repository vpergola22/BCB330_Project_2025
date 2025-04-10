# Load necessary libraries
library(immunarch)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(stringr)
library(tidyr)
library(gridExtra)


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



# ======= Clonal Proportion =======

imm_pr <- repClonality(immdata$data, .method = "clonal.prop", .clone.types = 0.5)
vis(imm_pr)


# ======= Overlap ========
# Compute and visualize repertoire overlap (adjusting public clonotype definition) --> general
repOverlap(immdata$data, .method = "public", .col = "v") %>% vis() + ggtitle("Repertoire Overlap of Public Clonotypes (v)") #TODO testing this out

repOverlap(immdata$data, .method = "public", .col = "j") %>% vis() + ggtitle("Repertoire Overlap of Public Clonotypes (j)")

repOverlap(immdata$data, .method = "public", .col = "v+j") %>% vis() + ggtitle("Repertoire Overlap of Public Clonotypes (v+j)")

repOverlap(immdata$data, .method = "public", .col = "aa+v+j") %>% vis() + ggtitle("Repertoire Overlap of Public Clonotypes (aa+v+j)")

repOverlap(immdata$data, .method = "public", .col = "aa") %>% vis() + ggtitle("Repertoire Overlap of Public Clonotypes (aa)")

repOverlap(immdata$data, .method = "public", .col = "nt") %>% vis() + ggtitle("Repertoire Overlap of Public Clonotypes (nt)")

#repOverlap(immdata$data, .method = "public", .col = "aa+j") %>% vis() + ggtitle("Repertoire Overlap of Public Clonotypes (aa+j)")


# ======= Distribution of Public Clonotype Frequencies =======
pubRep(immdata$data, .col = "v", .quant = "count", .coding = TRUE, .verbose = FALSE) %>% vis() + ggtitle("Distribution of Public Clonotype Frequencies (v)")

pubRep(immdata$data, .col = "j", .quant = "count", .coding = TRUE, .verbose = FALSE) %>% vis() + ggtitle("Distribution of Public Clonotype Frequencies (j)")

pubRep(immdata$data, .col = "aa+v", .quant = "count", .coding = TRUE, .verbose = FALSE) %>% vis() + ggtitle("Distribution of Public Clonotype Frequencies (aa+v)")

pubRep(immdata$data, .col = "aa+j", .quant = "count", .coding = TRUE, .verbose = FALSE) %>% vis() + ggtitle("Distribution of Public Clonotype Frequencies (aa+j)")

pubRep(immdata$data, .col = "v+j", .quant = "count", .coding = TRUE, .verbose = FALSE) %>% vis() + ggtitle("Distribution of Public Clonotype Frequencies (v+j)")

pubRep(immdata$data, .col = "aa", .quant = "count", .coding = TRUE, .verbose = FALSE) %>% vis() + ggtitle("Distribution of Public Clonotype Frequencies (aa)")



# ======= Number of Clonotypes ======= #TODO CLONES

# # Assuming that the patient identifier is stored in the names of the list items
# combined_data <- immdata$data %>%
#   purrr::imap_dfr(~ .x %>% 
#                     mutate(Patient = .y)) %>%  # Add Patient based on list name
#   mutate(raw_clonotype_id = strsplit(as.character(raw_clonotype_id), ";")) %>%
#   unnest_longer(raw_clonotype_id)

# # Check the structure of the combined data
# head(combined_data)

# # Count clonotypes per patient
# clonotype_counts <- combined_data %>%
#   group_by(Patient) %>%
#   summarise(Clonotype_Count = n_distinct(raw_clonotype_id))

# # Plot the count of unique clonotypes per patient
# ggplot(clonotype_counts, aes(x = Patient, y = Clonotype_Count, fill = Patient)) +
#   geom_bar(stat = "identity") +
#   labs(x = "Patient", y = "Number of Unique Clonotypes") +
#   theme_minimal()

# TODO attempting another method
tcell_counts <- repExplore(immdata$data, .method = "clones")
vis(tcell_counts)



# ======== Gene Usage ========

# Compute V gene usage for human TCR beta chains
#v_gene_usage <- geneUsage(immdata$data, "hs.trbv", .norm = T)

# Compute J gene usage for human TCR beta chains
#j_gene_usage <- geneUsage(immdata$data, "hs.trbj", .norm = T)

# Histogram of V gene usage
#vis(v_gene_usage, .by = "Sample", .meta = immdata$meta, .plot = "box")
# Boxplot of J gene usage
#vis(j_gene_usage, .by = "Sample", .meta = immdata$meta, .plot = "box")

# Perform PCA on V gene usage data
#v_gene_pca <- geneUsageAnalysis(v_gene_usage, .method = "js+pca+kmeans", .perp = .01, .verbose = F)

# Visualize PCA results
#vis(v_gene_pca, .plot = "clust")


# ======= Number of Clones =======
# Summing the number of clones per sample
total_tcells <- immdata$data %>%
  bind_rows(.id = "Sample") %>%  # Combine all samples into one dataframe with Sample ID
  group_by(Sample) %>%
  tally(name = "Total_T_Cells")  # Count the total number of clones per sample

# Print the result as a bar graph
ggplot(total_tcells, aes(x = Sample, y = Total_T_Cells, fill = Sample)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  labs(title = "Total T Cells Per Sample",
       x = "Sample",
       y = "Total T Cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels



# ======= Diversity =======
# Compute and visualize diversity grouped by sample TODO this looks so weird rn 
diversity_stats <- repDiversity(immdata$data, .col = "v+j")
vis(diversity_stats, .by = "Sample", .meta = immdata$meta)



# ======= Hyperexpanded =======
# Check if immdata$data exists and is not empty
if (!is.null(immdata$data) && length(immdata$data) > 0) {
  # Detect hyperexpanded clones and visualize
  map(names(immdata$data), ~ {
    dataset <- immdata$data[[.x]]  # Access each dataset by name
    
    if (!is.null(dataset) && nrow(dataset) > 0) {  # Ensure data is valid before applying repClonality
      repClonality(dataset, 
                   .method = "homeo",
                   .clone.types = c(Small = .001, Medium = .001, Large = .01, Hyperexpanded = 1)) %>% 
        vis() +
        ggtitle(paste("Hyperexpanded Clones - Patient:", .x))
    } else {
      print(paste("Warning: Empty or NULL dataset detected for", .x, "skipping..."))
      NULL
    }
  })
} else {
  print("Error: immdata$data is NULL or empty.")
}

  

# ======= Regression =======
# First, calculate percentage of hyperexpanded clones

# Compute hyperexpanded clones percentage for each sample
hyperexpanded_analysis <- map_dfr(names(immdata$data), ~ {
  dataset <- immdata$data[[.x]]
  
  # Compute total clones and hyperexpanded clones
  total_clones <- nrow(dataset)
  hyperexpanded_clones <- sum(dataset$Proportion > 0.01)  # Example threshold
  
  data.frame(
    Sample = .x,
    Total_Clones = total_clones,
    Hyperexpanded_Clones = hyperexpanded_clones,
    Percentage_Hyperexpanded = (hyperexpanded_clones / total_clones) * 100
  )
})

# Merge total T cells information
merged_data <- left_join(hyperexpanded_analysis, total_tcells, by = "Sample")

# Diagnostic print to check data
print("Merged Data Structure:")
print(str(merged_data))

# Regression model
if (nrow(merged_data) > 0) {
  # Fit regression model
  regression_model <- lm(Percentage_Hyperexpanded ~ Total_T_Cells, data = merged_data)
  
  # Summary of regression
  print("Regression Summary:")
  print(summary(regression_model))
  
  # Compute residuals
  merged_data$Residuals <- resid(regression_model)
  
  # Visualization of regression
  ggplot(merged_data, aes(x = Total_T_Cells, y = Percentage_Hyperexpanded)) +
    geom_point() +
    geom_smooth(method = "lm") +
    labs(
      title = "Regression: Total T Cells vs Percentage of Hyperexpanded Clones",
      x = "Total T Cells",
      y = "Percentage of Hyperexpanded Clones"
    ) +
    theme_minimal()
  
  # Optional: If you have additional variables like TMB
  # Uncomment and modify as needed
  # if (!is.null(merged_data$Tumor_Mutational_Burden)) {
  #   correlation_result <- cor(merged_data$Residuals, merged_data$Tumor_Mutational_Burden, method = "pearson")
  #   print("Correlation of Residuals with TMB:")
  #   print(correlation_result)
  # }
}


# ======= Num Shared Samples - Clonotype Frequency =======

library(immunarch)
library(dplyr)
library(ggplot2)
library(gridExtra)

visualize_clonotype_sharing <- function(immdata, methods = c("v", "j", "aa", "v+j")) {
  # Validate input
  if (is.null(immdata$data) || length(immdata$data) == 0) {
    stop("Invalid immunarch data: empty or missing data")
  }
  
  # Prepare a list to store visualization results
  visualization_plots <- list()
  
  # Loop through each method
  for (method in methods) {
    tryCatch({
      # Compute public repertoire with error handling
      pub_rep <- pubRep(immdata$data, .col = method, .quant = "count", .coding = TRUE, .verbose = TRUE)
      
      # Check if pub_rep is empty
      if (is.null(pub_rep) || nrow(pub_rep) == 0) {
        warning(paste("No data found for method:", method))
        next
      }
      
      # Convert to dataframe
      pub_rep_df <- as.data.frame(pub_rep)
      
      # Prepare data for analysis
      samples_per_clonotype <- pub_rep_df %>%
        # Count number of samples with non-zero values for each clonotype
        mutate(Num_Samples_Shared = rowSums(select(., starts_with("X")) > 0)) %>%
        # Compute total frequency across all samples
        mutate(Total_Frequency = rowSums(select(., starts_with("X"))))
      
      # Add rownames as a column if they exist
      if(!is.null(rownames(pub_rep_df))) {
        samples_per_clonotype$Clonotype <- rownames(pub_rep_df)
      }
      
      # Filter out rows with zero frequency
      samples_per_clonotype <- samples_per_clonotype %>%
        filter(Total_Frequency > 0)
      
      # Perform correlation analysis
      if (nrow(samples_per_clonotype) > 1) {
        correlation <- cor(samples_per_clonotype$Num_Samples_Shared, 
                           samples_per_clonotype$Total_Frequency, 
                           method = "spearman")
        cor_test <- cor.test(samples_per_clonotype$Num_Samples_Shared, 
                              samples_per_clonotype$Total_Frequency, 
                              method = "spearman")
        p_value <- cor_test$p.value
        
        # 1. Scatter plot with correlation
        scatter_plot <- ggplot(samples_per_clonotype, 
                               aes(x = Num_Samples_Shared, y = Total_Frequency)) +
          geom_point(alpha = 0.5, color = "darkblue") +
          geom_smooth(method = "lm", se = TRUE, color = "red") +
          labs(
            title = paste("Clonotype Sharing Analysis -", method),
            subtitle = sprintf("Spearman Correlation: %.3f (p-value: %.3f)", correlation, p_value),
            x = "Number of Samples Sharing Clonotype",
            y = "Total Clonotype Frequency"
          ) +
          theme_minimal()
        
        # 2. Box plot of frequency by number of shared samples
        box_plot <- ggplot(samples_per_clonotype, 
                           aes(x = factor(Num_Samples_Shared), y = Total_Frequency)) +
          geom_boxplot(fill = "lightblue") +
          labs(
            title = paste("Frequency Distribution by Sample Sharing -", method),
            x = "Number of Samples Shared",
            y = "Total Clonotype Frequency"
          ) +
          theme_minimal()
        
        # 3. Histogram of sample sharing
        histogram <- ggplot(samples_per_clonotype, aes(x = Num_Samples_Shared)) +
          geom_histogram(binwidth = 1, fill = "green", alpha = 0.7) +
          labs(
            title = paste("Distribution of Sample Sharing -", method),
            x = "Number of Samples Sharing Clonotype",
            y = "Count of Clonotypes"
          ) +
          theme_minimal()
        
        # Store plots
        visualization_plots[[method]] <- list(
          scatter_plot = scatter_plot,
          box_plot = box_plot,
          histogram = histogram,
          correlation = correlation,
          p_value = p_value,
          data = samples_per_clonotype
        )
      } else {
        warning(paste("Not enough data for correlation analysis in method:", method))
      }
    }, error = function(e) {
      warning(paste("Error processing method", method, ":", e$message))
    })
  }
  
  return(visualization_plots)
}

# Function to create and display plot grid
create_plot_grid <- function(visualizations) {
  if (length(visualizations) == 0) {
    stop("No visualizations to display")
  }
  
  # Create a list of plot grids
  plot_list <- lapply(names(visualizations), function(method) {
    method_plots <- visualizations[[method]]
    grid.arrange(
      method_plots$scatter_plot, 
      method_plots$box_plot, 
      method_plots$histogram, 
      ncol = 3,
      top = method
    )
  })
  
  return(plot_list)
}

# Example usage
# Ensure you load your immunarch data first
# immdata <- load_immunarch_data()

sharing_visualizations <- visualize_clonotype_sharing(immdata)
plot_grids <- create_plot_grid(sharing_visualizations)

# Generate and display the plot grid
plot_grids <- create_plot_grid(sharing_visualizations)