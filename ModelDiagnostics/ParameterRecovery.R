# Set working directory and source initialization file
setwd("~/GitDir/CodeWithPapers/PatternCompletion_CueFeatures/")
source("Initialise.R")
source("DDMFunctions/SimulateAccRT.R")
source("DDMFunctions/RunDDM.R")

############## Read in drift diffusion data and organise it ############

# Initialize empty dataframe to store all participants' data
driftdiffusion_data <- c()

# Loop through all participants and read their drift diffusion parameter files
for(part in unique(participants$Participant)){
  # Read TSV file containing drift diffusion parameters for each participant
  part_data <- read.table(paste0(data_path, part, "/", part, "_DriftDiffusionParameters.tsv"),
                          header=TRUE,
                          sep = "\t")
  # Combine with existing data
  driftdiffusion_data <- bind_rows(driftdiffusion_data, part_data)
}

# Prepare test data for analysis
# Set factor levels for proper statistical modeling
driftdiffusion_data <- driftdiffusion_data %>% 
  # Convert variables to factors with specified levels for proper statistical analysis
  mutate(Degradation = as.character(Degradation),
         Group = factor(Group,
                        levels = factor_labels$Group$levels),
         TrialType = factor(TrialType,
                            levels = factor_labels$TrialType$levels),
         Scrambling = factor(Scrambling,
                             levels = factor_labels$Scrambling$levels))

driftdiffusion_data_long <- driftdiffusion_data %>% 
  pivot_longer(cols = c(alpha, beta, tau, delta),
               names_to = "Parameter",
               values_to = "ParamVal")


ddparams_model8 <- bind_rows(
  driftdiffusion_data_long %>% 
    filter(Parameter %in% c("delta", "tau")),
  driftdiffusion_data_long %>% 
    filter(Parameter %in% c("alpha", "beta")) %>%
    group_by(Participant, Parameter, Group) %>% 
    summarise(ParamVal = mean(ParamVal)) %>% 
    mutate(Degradation = "collapsed",
           TrialType = "collapsed",
           Scrambling = "collapsed"))

############## Read in behavioural data (for comparison) and organise it ############

# Initialize empty dataframe to store all participants' data
behavioural_data <- c()

# Loop through all participants and read their task event files
for(part in unique(participants$Participant)){
  # Read TSV file containing task events for each participant
  part_data <- read.table(paste0(data_path, part, "/", part, "_TrialData.tsv"),
                          header=TRUE,
                          sep = "\t")
  # Combine with existing data
  behavioural_data <- bind_rows(behavioural_data, part_data)
}

# Filter to only include Test phase trials
test_data <- behavioural_data %>% 
  filter(Phase == "Test")

# Prepare test data for analysis
# Create accuracy variable and set factor levels for proper statistical modeling
test_data <- test_data %>% 
  # Determine what the correct response should be based on trial type
  # "Old" trials should be responded to as "Old", all others as "New"
  mutate(CorrectResponse = ifelse(TrialType == "Old",
                                  "Old", 
                                  "New"),
         # Calculate accuracy: compare participant's response to correct response
         Accuracy = CorrectResponse == Response,
         # Convert variables to factors with specified levels for proper statistical analysis
         Group = factor(Group,
                        levels = factor_labels$Group$levels),
         TrialType = factor(TrialType,
                            levels = factor_labels$TrialType$levels),
         Scrambling = factor(Scrambling,
                             levels = factor_labels$Scrambling$levels)) %>% 
  # Exclude bad trials (missing RT or RT <= 100ms) before simulating
  filter(!(is.na(RT) | RT <= 100)) %>% 
  # Remove the accuracy column. It will be recomputed following simulations
  select(-Accuracy)

############## Parameter Recovery ##############

# ---- Setup paths to save outputs ----
# Create output directory for parameter recovery results
param_recovery_path <- paste0(derivatives_path, "ParameterRecovery/")
if (!dir.exists(param_recovery_path)) {
  dir.create(param_recovery_path, recursive = TRUE)
}

# ---- Step 1: Simulate data using original DDM parameters ----
# For each participant, use their fitted DDM parameters (from Model 8)
# to generate simulated response data via the Wiener diffusion model.
# The simulation matches the number of trials per condition from the
# participant's actual empirical data, so the simulated dataset has the
# same structure as the original.
simulated_data <- lapply(
  unique(ddparams_model8$Participant),
  function(part) {
    cat(sprintf("Simulating data for participant: %s\n", part))

    # Get this participant's DDM parameters from Model 8.
    # For alpha/beta: one value per participant (collapsed across
    #   Degradation, Scrambling, TrialType).
    # For delta/tau: one value per condition (varies across
    #   Degradation, Scrambling, TrialType).
    part_params <- ddparams_model8 %>%
      filter(Participant == part)

    # Get this participant's empirical test data. This is used as the
    # comparison dataset to determine (a) which conditions to simulate
    # and (b) how many trials to simulate per condition.
    part_test_data <- test_data %>%
      filter(Participant == part)

    # Simulate responses for all conditions this participant saw
    sim_result <- simulate_DDM(
      mean_params_df = part_params,
      comparison_data = part_test_data
    )

    sim_result$simulations
  }
)

# Combine all participants' simulated data into a single dataframe
simulated_data <- do.call("rbind", simulated_data)

# Save the simulated trial-level data
write.csv(simulated_data,
          file = paste0(param_recovery_path, "SimulatedData.csv"),
          row.names = FALSE)

# ---- Step 2: Re-estimate DDM parameters from simulated data ----
# Keep only the simulated trials (discard the empirical copies that
# simulate_DDM returns alongside). These simulated trials have the
# same columns as the original test_data (Participant, Group,
# Degradation, Scrambling, TrialType, Response, RT).
simulated_trials <- simulated_data %>%
  filter(EmpSim == "Simulated")

# Fit DDM on the simulated data to attempt to recover the original
# parameters. If the model and fitting procedure are reliable, the
# recovered parameters should closely match ddparams_model8.
# Results are saved to derivatives_path/ParameterRecovery/.
recovered_results <- run_ddm(simulated_trials, "ParameterRecovery")

# Extract the recovered per-participant DDM parameters
recovered_params <- recovered_results$participant_ddm_params

# ---- Save simulation and recovery outputs ----
# Save the recovered DDM parameters
write.csv(recovered_params,
          file = paste0(param_recovery_path, "RecoveredParams.csv"),
          row.names = FALSE)

# ---- Step 3: Compare recovered vs original parameters ----
# The goal is to assess whether the DDM fitting procedure can reliably
# recover known parameter values. We correlate the original parameters
# (used to generate the simulated data) with the recovered parameters
# (re-estimated from the simulated data).

# Pivot recovered parameters from wide format (one column per parameter)
# to long format (Parameter + Recovered columns) to match ddparams_model8.
# Convert condition columns to character for consistent joining.
recovered_params_long <- recovered_params %>%
  pivot_longer(cols = c(alpha, beta, delta, tau),
               names_to = "Parameter",
               values_to = "Recovered")

# Match the structure of the original Model 8 parameters:
# In Model 8, alpha/beta are fixed across conditions (collapsed), while
# delta/tau vary per condition. To make a fair comparison, average the
# recovered alpha/beta across conditions before comparing.
recovered_for_comparison <- bind_rows(
  # delta and tau: keep per-condition values as-is
  recovered_params_long %>%
    filter(Parameter %in% c("delta", "tau")) %>% 
    mutate(Degradation = as.character(Degradation)),
  # alpha and beta: average across all conditions per participant,
  # since the original values were also collapsed
  recovered_params_long %>%
    filter(Parameter %in% c("alpha", "beta")) %>%
    group_by(Participant, Parameter, Group) %>%
    summarise(Recovered = mean(Recovered), .groups = "drop") %>%
    mutate(Degradation = "collapsed",
           TrialType = "collapsed",
           Scrambling = "collapsed")
)

# Prepare original parameters for joining by ensuring consistent types
original_for_comparison <- ddparams_model8 %>%
  rename(Original = ParamVal)

# Join original and recovered parameters on all matching keys
param_comparison <- full_join(
  original_for_comparison,
  recovered_for_comparison,
  by = c("Participant", "Parameter", "Group",
         "Degradation", "Scrambling", "TrialType")
)

# Compute Pearson's correlation between original and recovered values,
# separately for each DDM parameter
param_correlations <- param_comparison %>%
  group_by(Parameter) %>%
  cor_test(Original, Recovered, method = "pearson")

print(param_correlations)

# Scatterplot of original vs recovered parameters, faceted by parameter.
# The dashed identity line (y = x) indicates perfect recovery.

xy_range <- c(min(min(param_comparison$Recovered, na.rm = T), 
                  min(param_comparison$Original, na.rm = T)),
              max(max(param_comparison$Recovered, na.rm = T), 
                  max(param_comparison$Original, na.rm = T)))
recovery_scatter <- ggscatter(
  param_comparison,
  x = "Original",
  y = "Recovered",
  size = 3,
  alpha = 0.5,
  color = "Parameter",
) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed") +
  labs(x = "Original",
       y = "Recovered") +
  coord_cartesian(xlim = xy_range,
                  ylim = xy_range) +
  x_axis_theme + y_axis_theme +
  blank_bg_theme + paper_facet_theme

# Save the scatterplot and correlation results
png(paste0(param_recovery_path, "ParameterRecovery_Scatter.png"),
    width = 1200, height = 1000)
print(recovery_scatter)
dev.off()

# Scatterplot of original vs recovered parameters, coloured by parameter.
# Only plotting first 6 participants for visualisation
# The dashed identity line (y = x) indicates perfect recovery.

parts_to_plot <- participants$Participant[1:6]
param_comparison_parts <- param_comparison %>% 
  filter(Participant %in% parts_to_plot)

xy_range <- c(min(min(parts_to_print$Recovered, na.rm = T), 
                  min(parts_to_print$Original, na.rm = T)),
              max(max(parts_to_print$Recovered, na.rm = T), 
                  max(parts_to_print$Original, na.rm = T)))
recovery_scatter <- ggscatter(
  param_comparison,
  x = "Original",
  y = "Recovered",
  size = 3,
  alpha = 0.5,
  color = "Parameter",
) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed") +
  labs(x = "Original",
       y = "Recovered") +
  coord_cartesian(xlim = xy_range,
                  ylim = xy_range) +
  x_axis_theme + y_axis_theme +
  blank_bg_theme + paper_facet_theme

write.csv(param_correlations,
          file = paste0(param_recovery_path,
                        "ParameterCorrelations.csv"),
          row.names = FALSE)
