# Load required libraries
if (!requireNamespace("immunarch", quietly = TRUE)) install.packages("immunarch")
library(immunarch)
library(dplyr)
library(readr)
library(stringr)
library(SingleCellExperiment)
library(tidySingleCellExperiment)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(pheatmap)
library(ggraph)
library(tidygraph)

# Define base directories
data(immdata)
  
  
  
    
    # Now proceed with the analysis
    names(immdata)
    names(immdata$meta)
    
    # doesnt work !!!! same error i was getting on my data set as well
    # basic_stats <- repExplore(immdata, .method = "volume")
    # basic_stats



    #immdata$meta
    #imm_pr <- repClonality(immdata$data, .method = "clonal.prop")
    #imm_pr
    # write_csv(exp, file.path("/ddn_exa/campbell/vpergola/Data/Pipeline/immunarch/song2022/basic_statistics.csv"))

