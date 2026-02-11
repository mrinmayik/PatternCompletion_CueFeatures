# Set working directory and source initialization file
setwd("~/GitDir/CodeWithPapers/PatternCompletion_CueFeatures/")
source("Initialise.R")
source("RunDDM.R")

############## Read in behavioural data and organise it ############

# Initialize empty dataframe to store all participants' data
behavioural_data <- c()

# Loop through all participants and read their task event files
for(part in unique(participants$Participant)){
  # Read TSV file containing task events for each participant
  part_data <- read.table(paste0(data_path, part, "/", part, "_Data.tsv"),
                          header = TRUE,
                          sep = "\t")
  # Combine with existing data
  behavioural_data <- bind_rows(behavioural_data, part_data)
}

# Filter to only include Test phase trials
test_data <- behavioural_data %>% 
  filter(Phase == "Test")

# Fit the diffusion model by calling the run_ddm function
run_success <- run_ddm(test_data, "MaximalModel")
