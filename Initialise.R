#!/bin/bash Rscript

################## Package Management ##################
# Check if required packages are installed and install if missing
required_packages <- c("stringr", "tidyr", "dplyr", "ggplot2", "ggpubr", "ez", "rstatix",
                       "psychReport")
to_be_installed <- required_packages[!(required_packages %in% rownames(installed.packages()))]
if(length(to_be_installed) > 0){
  for(package in to_be_installed){
    install.packages(package)
  }
}

################## Library Imports ##################
# Load all required libraries for data analysis and visualization
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ez)
library(rstatix)
library(psychReport)

################## Path Setup [UPDATE THESE FOR YOUR LOCAL MACHINE] ##################
# Define base paths for scripts and data
scripts_path <- "~/GitDir/CodeWithPapers/PatternCompletion_CueFeatures/"

# Change base_path to the folder that contains the downloaded data from OpenNeuro
# All other paths are defined relative to base_path
base_path <- "~/Desktop/PatternCompletion/"
data_path <- paste0(base_path, "DataForSharing/")
derivatives_path <- paste0(data_path, "derivatives/")

################## Participant Management ##################
# Read and filter participant information based on task completion
participants <- read.table(paste0(data_path,
                                  "/participants.tsv"),
                           header = TRUE)

################## Factor Level Definitions ##################
# Define factor levels and labels for different experimental conditions
# Each list element contains standardized levels, display labels, and (where applicable) color codes

factor_labels <- list(
  "Degradation" = 
    list(levels = c(40, 70, 85, 95),
         axis_label = "Degradation"),
  "Scrambling" = 
    list(levels = c("Intact", "Scrambled"),
         axis_label = "Scrambling"),
  "TrialType" = 
    list(levels = c("Old", "Similar", "New"),
         axis_label = "Trial Type"),
  "Group" = 
    list(levels = c("YA", "OA"),
         axis_label = "Group",
         colours = c("maroon", "steelblue"))
)

################## Utility Functions ##################
# Function to correct p-values using different methods
correct_p_vals <- function(df, pvalcol, usemethod = "fdr"){
  df <- as.data.frame(df)
  #Is raw p-value significant?
  df[, "sig"] <- ifelse(df[, pvalcol] <= 0.05, "*", "")
  
  #Correct p's based on method entered
  df[, paste(pvalcol, "_", usemethod, "corrected", sep = "")] <- p.adjust(df[, pvalcol], method = usemethod)
  #Is corrected p significant?
  df[, paste("sig_", usemethod, "corrected", sep = "")] <- ifelse(df[, paste(pvalcol, "_", usemethod, "corrected", sep = "")] <= 0.05, "*", "")
  
  return(df)
}

# Function to calculate summary statistics for a given column
summarise_data <- function(df, col_name, rm_na = FALSE){
  df <- as.data.frame(df)
  M <- mean(df[,col_name], na.rm = rm_na)
  SD <- sd(df[,col_name], na.rm = rm_na)
  SE <- SD / sqrt(nrow(df))
  LCI <- M - 1.96*SE
  HCI <- M + 1.96*SE
  MeanPlusSE <- M+SE
  MeanMinusSE <- M-SE
  NumOfRows <- nrow(df)
  data.frame(Mean=M, SD=SD, SE=SE, LCI=LCI, HCI=HCI, MeanPlusSE=MeanPlusSE, MeanMinusSE, Rows=NumOfRows)
}

################## Plotting Themes ##################
# Define consistent themes for paper-quality visualizations
# Theme for facet labels and panel borders
paper_facet_theme <- theme(strip.text.x = element_text(size = 22, 
                                                       colour = "black"),
                           strip.text.y = element_text(size = 22, 
                                                       colour = "black"), 
                           strip.background = element_rect(color = "white", 
                                                           fill = "white", 
                                                           linewidth = 1.5, 
                                                           linetype = "solid"),
                           panel.border = element_rect(colour = "black", 
                                                       fill = NA, 
                                                       linewidth=1.5),
                           strip.placement = "outside",)
# Theme for x-axis formatting
x_axis_theme <- theme(axis.title.x = element_blank(), 
                      axis.text.x = element_text(colour = "#000000",
                                                 size = 20)) 
# Theme for y-axis formatting
y_axis_theme <-   theme(axis.title.y = element_text(colour = "#000000",
                                                    size = 22), 
                        axis.text.y = element_text(colour = "#000000",
                                                   size = 20))
# Theme for clean background without grid lines
blank_bg_theme <- theme(panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(), 
                        axis.line = element_line(colour = "black"))
# Theme for legend formatting
legend_theme <- theme(legend.text = element_text(face = "bold", size = 20), 
                      legend.title = element_text(face = "bold", size = 25),
                      legend.position = "bottom")
# Theme for title formatting
title_theme <- theme(plot.title = element_text(face = "bold", 
                                               size = 25, 
                                               hjust = 0.5), 
                     legend.key.size=unit(1.3, "cm"))


